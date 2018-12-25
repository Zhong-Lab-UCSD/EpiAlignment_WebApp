const ClusterProcesser = require('./clusterProcessing')

module.exports = function () {
  this.clusterProc = new ClusterProcesser(
    [
      {
        name: 'human',
        latin: 'Homo_sapiens',
        reference: 'hg38',
        encode_reference: 'GRCh38'
      }, {
        name: 'mouse',
        latin: 'Mus_musculus',
        reference: 'mm10',
        encode_reference: 'mm10'
      }
    ],
    {
      rawFilePath: 'epialignment_beta/Annotation/AnnotationFiles',
      clusterSuffix: '_clusters'
    }
  )
}
