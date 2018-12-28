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

const geneAnnoPath = 'Annotation/AnnotationFiles/genes'
const rawClusterFilePath = 'Annotation/AnnotationFiles'
const clusterSuffix = '.gene_info.gz'

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24
const MAX_MATCH_ENTRIES = 10

const serverBasePath = '.'

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
  constructor (symbol, ensemblId, aliases, description) {
    this.symbol = symbol
    this.ensemblId = ensemblId
    this.aliases = aliases || []
    this.aliases.unshift(this.symbol.toLowerCase())
    this.description = description || ''
  }

  toJSON () {
    return {
      symbol: this.symbol,
      aliases: this.aliases,
      ensemblId: this.ensemblId,
      description: this.description
    }
  }

  static getGeneFromNcbiEntry (ncbiEntry) {
    let tokens = ncbiEntry.trim().split('\t')
    let ensemblId = null
    tokens[5].split('|') // dbXrefs, use to populate Ensembl ID
      .some(entry => {
        let [key, value] = entry.split(/:(.+)/)
        if (key === 'Ensembl') {
          ensemblId = value
        }
      })
    let description = tokens[8]
    let symbol = tokens[2]
    let aliases = tokens[4] !== '-' ? tokens[4].split('|') : []
    return new this(symbol, ensemblId, aliases, description)
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
  constructor (id, geneEnsemblId, speciesName, geneMap, geneName) {
    this.id = id
    this.genesBySpecies = {}
    this.addNewSpeciesIfNotExists(speciesName)
    this.addGeneById(geneEnsemblId, speciesName, geneMap, geneName)
  }

  addNewSpeciesIfNotExists (speciesName) {
    if (!this.genesBySpecies.hasOwnProperty(speciesName)) {
      this.genesBySpecies[speciesName] = []
      this.genesBySpecies[speciesName]._map = new Map()
    }
  }

  addGeneById (geneEnsemblId, speciesName, geneMap, geneName) {
    this.addNewSpeciesIfNotExists(speciesName)
    if (!geneMap[speciesName].has(geneEnsemblId)) {
      // No such entry from NCBI data base, build a temp one
      geneMap[speciesName].set(geneEnsemblId,
        new Gene(geneName, geneEnsemblId))
    }
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
  let data = String(await zlipGunzipPromise(buffer))
  let speciesGeneMap = new Map()
  data.split('\n').forEach(line => {
    if (line && !line.startsWith('#')) {
      let newGene = Gene.getGeneFromNcbiEntry(line)
      if (newGene.ensemblId) {
        speciesGeneMap.set(newGene.ensemblId, newGene)
      }
    }
  })
  return speciesGeneMap
}

function getSpeciesGeneAnnoFileName (species) {
  return path.format({
    dir: path.format({ dir: serverBasePath, base: geneAnnoPath }),
    name: species.latin,
    ext: geneInfoSuffix
  })
}

async function loadSpeciesGeneAnnoFromUrlOrBuffer (species) {
  let geneAnnoFileName = getSpeciesGeneAnnoFileName(species)
  let needsUpdate = await fsStatPromise(geneAnnoFileName)
    .catch(err => {
      if (err.code === 'ENOENT') {
        return true
      } else {
        console.log(err)
        throw err
      }
    }).then(geneAnnoFileStat =>
      (Date.now() - geneAnnoFileStat.ctime > (90 * MILLISECONDS_IN_A_DAY)))
  let buffer
  if (needsUpdate) {
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
    dir: settings.rawFilePath || path.format({
      dir: serverBasePath,
      base: rawClusterFilePath
    }),
    name: species.reference,
    ext: settings.clusterSuffix || clusterSuffix
  })
  let buffer = await fsReadFilePromise(clusterFileName, 'utf8')
  buffer.split('\n').forEach(line => {
    if (line) {
      let clusterId = Cluster.getIdFromLine(line)
      let geneEnsemblId = line.split('\t')[0]
      let geneName = line.split('\t')[1]
      // console.log(geneEnsemblId)
      if (!clusters._map.has(clusterId)) {
        let newCluster =
          new Cluster(
            clusterId, geneEnsemblId, species.name, geneMap, geneName)
        clusters.push(newCluster)
        clusters._map.set(clusterId, newCluster)
      } else {
        clusters._map.get(clusterId)
          .addGeneById(geneEnsemblId, species.name, geneMap, geneName)
      }
    }
  })
  return clusters
}

function addAliasClusterIdPair (alias, clusterId, mapObj, keepCase) {
  alias = alias || ''
  if (!keepCase) {
    alias = alias.toLowerCase()
  }
  if (!mapObj.has(alias)) {
    mapObj.set(alias, clusterId)
  } else {
    if (Array.isArray(mapObj.get(alias))) {
      if (
        mapObj.get(alias).indexOf(clusterId) < 0
      ) {
        mapObj.get(alias).push(clusterId)
      }
    } else {
      if (mapObj.get(alias) !== clusterId) {
        mapObj.set(alias, [
          mapObj.get(alias), clusterId
        ])
      }
    }
  }
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
    this.ensemblIdToClusterMap = new Map()
    console.log('===== Building cluster table =====')
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
            addAliasClusterIdPair(
              gene.ensemblId, cluster.id, this.ensemblIdToClusterMap, true
            )
            gene.aliases.forEach(alias =>
              addAliasClusterIdPair(
                alias, cluster.id, this.aliasToClusterMap
              )
            )
          })
        }
      }
    })).then(() => console.log('===== Cluster table built ====='))
  }

  getClusters (partialAlias, maxMatchEntries) {
    maxMatchEntries = maxMatchEntries || MAX_MATCH_ENTRIES
    if (partialAlias.startsWith('ENSG') ||
      partialAlias.startsWith('ENSMUSG') ||
      partialAlias.match(/ENS[A-Z]+[0-9]{3}/)
    ) {
      // ensembl ID
      return this.clusterTableReadyPromise.then(() => {
        let fullMatchList = this.ensemblIdToClusterMap.has(partialAlias)
          ? this.ensemblIdToClusterMap.get(partialAlias)
          : []
        if (!Array.isArray(fullMatchList)) {
          fullMatchList = [fullMatchList]
        }
        return {
          fullMatchList: fullMatchList.map(
            clusterId => this.clusters._map.get(clusterId)),
          maxExceeded: false,
          partialMatchList: []
        }
      })
    } else {
      partialAlias = (partialAlias || '').toLowerCase()
      return this.clusterTableReadyPromise.then(() => {
        let fullMatchSet = new Set()
        let fullMatchList = this.aliasToClusterMap.has(partialAlias)
          ? this.aliasToClusterMap.get(partialAlias)
          : []
        if (!Array.isArray(fullMatchList)) {
          fullMatchList = [fullMatchList]
        }
        fullMatchList.forEach(cluster => fullMatchSet.add(cluster))
        let partialMatchMap = new Map()
        let maxExceeded = false
        for (const [key, value] of this.aliasToClusterMap) {
          if (key !== partialAlias && key.includes(partialAlias)) {
            let valueArray = Array.isArray(value) ? value : [value]
            valueArray.forEach(clusterId => {
              if (!fullMatchSet.has(clusterId) &&
                !partialMatchMap.has(clusterId)
              ) {
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
          partialMatchList: maxExceeded ? [] : [...partialMatchMap.values()]
        }
      }).catch(err => {
        console.log(err)
        return null
      })
    }
  }
}

module.exports = ClusterProcesser
