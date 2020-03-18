/* This is a small config file that specifies genomes
  and associated bait / VCF files */

//Select genome build
def genome="GRCh38" // Choose GRCh37d5 or GRCh38

//Select bait bait_capture_files ensure the corresponding build is specified
def bait_


// Scripts specifying genome singularity paths
if(genome == "GRCh38") {
  // defining genome path
  def core="hg38/UCSC/WholeGenome_hg38.tar.gz"
  def
  println "GRCh38 genome selected"
} else if(genome == "GRCh37d5") {
  def core="hg38/UCSC/WholeGenome_hg38.tar.gz"
  println "GRCh37d5 genome selected"
} else {
  println "Please choose correct genome."
}


// Scripts specifying bait paths
