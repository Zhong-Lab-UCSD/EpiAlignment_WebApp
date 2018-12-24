const axios = require('axios')
const path = require('path')
const fs = require('fs')
const util = require('util')
const zlib = require('zlib')

const fsReadFilePromise = util.promisify(fs.readFile)
const fsWriteFilePromise = util.promisify(fs.writeFile)
const fsStatPromise = util.promisify(fs.stat)
const zlipGunzipPromise = util.promisify(zlib.gunzip)

const ncbiUrl = 'https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/'
const geneInfoSuffix = '.gene_info.gz'

const geneAnnoPath = 'Annotation/AnnotationFiles/gene'
const rawClusterFilePath = 'Annotation/AnnotationFiles'
const clusterSuffix = '.gene_info.gz'

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24
const MAX_MATCH_ENTRIES = 10

/**
 * Workflow:
 *
 * When initializing
 * *  Build an NCBI query table for `Gene` objects, keyed by their Ensembl IDs
 *    (Check if the NCBI gene_info file needs to be updated.)
 * *  Read cluster file from `settings.rawFilePath` + `reference` +
 *    `settings.clusterPostfix`
 * *  Build a list of clusters, with all its genes inserted
 * *  Enum through the clusters and populate alias/ensemblId/symbol-to-gene
 *    map
 * *  Mark completion of initialization
 *
 * When querying
 * *  Wait until initialization is complete (or promise rejected)
 * *  Enum all possible names in the map and find everything that's partially
 *    matched or completely matched (in two categories)
 * *  Return the JSON for two categories of clusters
 */

/**
 * Gene object
 *
 * @class Gene
 * @property {string} symbol Gene symbol ("official gene name") as is defined
 *    by NCBI
 * @property {Array<string>} aliases Gene aliases, from NCBI (note that this
 *    includes `this.symbol` as well, unlike NCBI files)
 * @property {string} ensemblId Ensembl ID
 * @property {string} description Gene description
 *
 * @constructor
 * @param {string} ncbiEntry the entry in NCBI gene_info file
 */
class Gene {
  constructor (ncbiEntry) {
    let tokens = ncbiEntry.trim().split('\t')
    this.ensemblId = null
    tokens[5].split('|') // dbXrefs, use to populate Ensembl ID
      .some(entry => {
        let [key, value] = entry.split(/:(.+)/)
        if (key === 'Ensembl') {
          this.ensemblId = value
        }
      })
    this.description = tokens[7]
    this.symbol = tokens[2]
    this.aliases = tokens[4] !== '-'
      ? tokens[4].split('|').unshift(this.symbol)
      : [this.symbol]
  }
  toJSON () {
    return {
      symbol: this.symbol,
      aliases: this.aliases,
      ensemblId: this.ensemblId,
      description: this.description
    }
  }
}

/**
 * Class for a gene cluster
 *
 * @class Cluster
 * @property {object} genesBySpecies an object with species name as a key,
 *    the value is an Array of all genes included in this cluster in a given
 *    species.
 * @property {string} id Unique cluster ID
 */
class Cluster {
  constructor (id, geneEnsemblId, speciesName, geneMap) {
    this.id = id
    this.genesBySpecies = {}
    this.addNewSpeciesIfNotExists(speciesName)
    this.addGeneById(geneEnsemblId, speciesName, geneMap)
  }

  addNewSpeciesIfNotExists (speciesName) {
    if (!this.genesBySpecies.hasOwnProperty(speciesName)) {
      this.genesBySpecies[speciesName] = []
      this.genesBySpecies[speciesName]._map = new Map()
    }
  }

  addGeneById (geneEnsemblId, speciesName, geneMap) {
    this.addNewSpeciesIfNotExists(speciesName)
    if (!this.genesBySpecies[speciesName]._map.has(geneEnsemblId)) {
      this.genesBySpecies[speciesName].push(
        geneMap[speciesName].get(geneEnsemblId)
      )
      this.genesBySpecies[speciesName]._map.set(
        geneEnsemblId, geneMap[speciesName].get(geneEnsemblId))
    }
  }

  static getIdFromLine (line) {
    return line.split('\t')[2]
  }
}

async function loadGeneAnnoFromGzipBuffer (buffer) {
  let data = await zlipGunzipPromise(buffer)
  let speciesGeneMap = new Map()
  data.split('\n').forEach(line => {
    if (line && !line.startsWith('#')) {
      let newGene = new Gene(line)
      if (newGene.ensemblId) {
        speciesGeneMap.set(newGene.ensemblId, newGene)
      }
    }
  })
  return speciesGeneMap
}

function getSpeciesGeneAnnoFileName (species) {
  return path.format({
    dir: geneAnnoPath,
    name: species.latin,
    ext: geneInfoSuffix
  })
}

async function loadSpeciesGeneAnnoFromUrlOrBuffer (species) {
  let geneAnnoFileName = getSpeciesGeneAnnoFileName(species)
  let geneAnnoFileStat = await fsStatPromise(geneAnnoFileName)
  let buffer
  if (Date.now() - geneAnnoFileStat.ctime > (90 * MILLISECONDS_IN_A_DAY)) {
    // file is more than 90 days old, update from NCBI server
    buffer = (await axios.request({
      url: ncbiUrl + species.latin + geneInfoSuffix,
      responseType: 'arraybuffer'
    })).data
    fsWriteFilePromise(geneAnnoFileName, buffer)
  } else {
    buffer = await fsReadFilePromise(geneAnnoFileName)
  }
  return loadGeneAnnoFromGzipBuffer(buffer)
}

async function loadClustersFromFile (clusters, species, settings, geneMap) {
  let clusterFileName = path.format({
    dir: settings.rawFilePath || rawClusterFilePath,
    name: species.name,
    ext: settings.clusterSuffix || clusterSuffix
  })
  let buffer = await fsReadFilePromise(clusterFileName, 'utf8')
  buffer.split('\n').forEach(line => {
    if (line) {
      let clusterId = Cluster.getIdFromLine(line)
      let geneEnsemblId = line.split('\t')[0]
      if (!clusters._map.has(clusterId)) {
        let newCluster =
          new Cluster(clusterId, geneEnsemblId, species.name, geneMap)
        clusters.push(newCluster)
        clusters._map.set(clusterId, newCluster)
      } else {
        clusters._map.get(clusterId)
          .addGeneById(geneEnsemblId, species.name, geneMap)
      }
    }
  })
  return clusters
}

/**
 * A cluster processor for cluster search
 *
 * @class ClusterProcesser
 * @property {Array<Cluster>} clusters
 */
class ClusterProcesser {
  constructor (supportedSpecies, settings) {
    this.geneMap = {}
    this.clusters = []
    this.clusters._map = new Map()
    this.aliasToClusterMap = new Map()
    this.clusterTableReadyPromise = Promise.all(supportedSpecies.map(
      species => loadSpeciesGeneAnnoFromUrlOrBuffer(species)
        .then(speciesGeneMap => (this.geneMap[species.name] = speciesGeneMap))
    )).then(() => Promise.all(supportedSpecies.map(species =>
      loadClustersFromFile(this.clusters, species, settings, this.geneMap)
    ))).then(() => this.clusters.forEach(cluster => {
      for (let species in cluster.genesBySpecies) {
        if (cluster.genesBySpecies.hasOwnProperty(species)) {
          let geneList = cluster.genesBySpecies[species]
          geneList.forEach(gene => {
            gene.aliases.forEach(alias =>
              this.addAliasClusterIdPair(alias, cluster.id))
          })
        }
      }
    }))
  }

  addAliasClusterIdPair (alias, clusterId) {
    alias = (alias || '').toLowerCase()
    if (!this.aliasToClusterMap.has(alias)) {
      this.aliasToClusterMap.set(alias, clusterId)
    } else {
      if (Array.isArray(this.aliasToClusterMap.get(alias))) {
        if (
          this.aliasToClusterMap.get(alias).indexOf(clusterId) < 0
        ) {
          this.aliasToClusterMap.get(alias).push(clusterId)
        }
      } else {
        if (this.aliasToClusterMap.get(alias) !== clusterId) {
          this.aliasToClusterMap.set(alias, [
            this.aliasToClusterMap.get(alias), clusterId
          ])
        }
      }
    }
  }

  getClusters (partialAlias, maxMatchEntries) {
    maxMatchEntries = maxMatchEntries || MAX_MATCH_ENTRIES
    partialAlias = (partialAlias || '').toLowerCase()
    return this.clusterTableReadyPromise.then(() => {
      let fullMatchList = this.aliasToClusterMap.has(partialAlias)
        ? this.aliasToClusterMap.get(partialAlias)
        : []
      if (!Array.isArray(fullMatchList)) {
        fullMatchList = [fullMatchList]
      }
      let partialMatchMap = new Map()
      let maxExceeded = false
      for (const [key, value] of this.aliasToClusterMap) {
        if (key.includes(partialAlias)) {
          let valueArray = Array.isArray(value) ? value : [value]
          valueArray.forEach(clusterId => {
            if (!partialMatchMap.has(clusterId)) {
              partialMatchMap.set(clusterId,
                this.clusters._map.get(clusterId))
            }
          })
          if (partialMatchMap.size > maxMatchEntries) {
            maxExceeded = true
            break
          }
        }
      }
      return {
        fullMatchList: fullMatchList.map(
          clusterId => this.clusters._map.get(clusterId)),
        maxExceeded: maxExceeded,
        partialMatchList: maxExceeded ? [] : partialMatchMap.values()
      }
    }).catch(err => {
      console.log(err)
      return null
    })
  }
}

module.exports = ClusterProcesser
