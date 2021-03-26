function initBrowser() {

  var div,
          options,
          browser;

  div = document.getElementById("myDiv");
  options = {
      genome: {
    "id":"canFam2",
    "name": "canFam2",
    "fastaURL": "static/data/public/canFam2.fa",
    "indexURL": "static/data/public/canFam2.fa.fai",
    
    "tracks": [
      {
        "name": "Refseq Genes",
        "url": "https://s3.amazonaws.com/igv.org.genomes/hg38/refGene.txt.gz",
        "order": 1000000,
        "indexed": false
      }
    ]
  },
      //locus: "22:24,375,771-24,376,878",
      tracks: [
                  {
                    type: "variant",
                    format: "vcf",
                    url: "static/data/public/test.vcf.gz",
                    indexURL: "static/data/public/test.vcf.gz.tbi",
                    name: "test.vcf",
                    squishedCallHeight: 1,
                    expandedCallHeight: 10,
                    displayMode: "expanded",
                    visibilityWindow: 1000000000
                    }
              ]
  };

  browser = igv.createBrowser(div, options);
}
