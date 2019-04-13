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

const MILLISECONDS_IN_A_DAY = 1000 * 60 * 60 * 24

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

async function loadGeneAnnoFromGzipBuffer (buffer, keys, caseInSensitive) {
  keys = keys || ['ensemblId']
  if (!Array.isArray(keys)) {
    keys = [keys]
  }
  let data = String(await zlipGunzipPromise(buffer))
  let speciesGeneMap = new Map()
  data.split('\n').forEach(line => {
    if (line && !line.startsWith('#')) {
      let newGene = Gene.getGeneFromNcbiEntry(line)
      keys.forEach(key => {
        if (newGene.hasOwnProperty(key) && newGene[key]) {
          speciesGeneMap.set(
            caseInSensitive ? newGene[key].toLowerCase() : newGene[key],
            newGene)
        }
      })
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

async function loadSpeciesGeneAnno (species, keys, caseInSensitive) {
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
  return loadGeneAnnoFromGzipBuffer(buffer, keys, caseInSensitive)
}

module.exports.Gene = Gene
module.exports.loadSpeciesGeneAnno = loadSpeciesGeneAnno
