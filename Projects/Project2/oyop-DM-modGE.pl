#Exploring Bioinformatics: A Project-Based Approach
#Caroline St.Clair and Jonathan Visick

#Chapter 7 - On Your Own
# This implementation is for Neighbor Joining and only reads
# in distance matrices

# Modified significantly by Gunes 
# to allow distance matrices to be entered
# rather than distances to be computed!

%hashOriginal = (); # hash of hash to hold initial distance matrix value

%hashClusters = (); # hash of hash to hold matrix value
                    #    key = code for sequence 1
                    #    value = hash of sequence 2 code and distance

%hashNewick = (); # hash table to hold Newick code of cluster

my($data); # temp variable to hold current data
my(@codes); # This is the character code for each original OTUs
my($numOTU); # This is the number of OTUs (or the length of the matrix)
my(@distmatrix); # This is the "matrix" version of the distance matrix!
my(@temp); # temporary array

my($MAXDIST,$smallest,$smallestI,$smallestJ);

# Note that the format assumed for the distance matrix is similar to the 
# input substitution matrix format of previous projects
print "Enter the name of the distance matrix file: ";
$fileinput = <STDIN>;
chomp $fileinput;
if(!open(infile1, $fileinput)){
    print "error opening input file 1\n";
    exit;
}

# read the number of OTUs on line 1
$data = <infile1>;
chomp $data;
$numOTU = $data;

# now get the keys/codes of OTUs on line 2
$data = <infile1>;
chomp $data;
@codes = split(/[,\s+]/, $data); # this is how the keys are read

# Initialize hash codes and populate distance matrix
for ($i=0; $i<$numOTU; $i++) {
    $hashNewick{$codes[$i]} = $codes[$i];
    $hashOriginal{$codes[$i]} = $codes[$i];   
    $data = <infile1>;
    chomp $data;
    @temp = split(/[,\s+]/, $data);
    push @distmatrix, [@temp];
}

# Modified to Print out the distance matrix, initializing the hashClusters in the meantime
print "\n    ";
print "@codes\n";
for ($i = 0; $i < $numOTU; $i++) {
    print "$codes[$i]: ";
    for ($j = 0; $j < $numOTU; $j++) {
	print "$distmatrix[$i][$j] "; 
	$hashClusters{$codes[$i]}{$codes[$j]} = $distmatrix[$i][$j];
    }
    print "\n";
}

#########################################################
# Agglomerative Clustering Algorithm Implementing Neighbor-Joining Method
#     - merge clusters with smallest TRANSITION distance
#     - repeat process until two clusters remain
#########################################################
#@arrayClusters = keys(%hashOriginal); # initializing cluster array
$MAXDIST = 999999;
$numClusters = $numOTU; # initializing the number of clusters
while ($numClusters > 2){
  @arrayClusters = keys(%hashClusters);  
   # Calculate r values
   @arrayRValues = ();
   for ($i = 0; $i < $numClusters; $i++){
       $temp = 0;
       for ($j = 0; $j < $numClusters; $j++){
           #since value can be stored as hash{x}{y} or hash{y}{x}
           $temp = $temp
                     + $hashClusters{$arrayClusters[$i]}{$arrayClusters[$j]};
        }
        $arrayRValues[$i] = $temp / ($numClusters - 2); # with divide by 2 modification due to distances being calculated twice
   }
   print "R values\n";
   for ($i = 0; $i < $numClusters; $i++){
       print $arrayClusters[$i], " " , $arrayRValues[$i], " ";
   }
   print "\n";
   
   # Determine merging clusters using transformed distances.
  $smallest = $MAXDIST; $smallestI = 0; $smallestJ = 0;
  print "TD matrix: \n";
  for ($i = 0; $i < $numClusters-1; $i++){
      for ($j = $i+1; $j < $numClusters; $j++){
          # calculate transformed distance
          $tempTD = $hashClusters{$arrayClusters[$i]}{$arrayClusters[$j]}
                    - $arrayRValues[$i] - $arrayRValues[$j];
          # check if smallest transformed distance
	  print "TD($arrayClusters[$i],$arrayClusters[$j]) = $tempTD,   ";
          if ($tempTD < $smallest){
             $smallest = $tempTD;
             $smallestI = $i;
             $smallestJ = $j;
          }
      }
      print "\n"; 
   }
   # merge smallestI and smallestJ into clusters and add
   # to hash tables
   $clusterI = $arrayClusters[$smallestI];
   $clusterJ = $arrayClusters[$smallestJ];
   $merge = $clusterI.$clusterJ;
   # calculate branch lengths
   $branch1 = ($hashClusters{$clusterJ}{$clusterI}
              + $arrayRValues[$smallestI] - $arrayRValues[$smallestJ]) / 2;
   $branch2 = ($hashClusters{$clusterJ}{$clusterI}
              + $arrayRValues[$smallestJ] - $arrayRValues[$smallestI]) / 2;

   print "Merging Clusters: $clusterI and $clusterJ\n";
   print "    Distance between $clusterI and ancestral node = $branch1\n";
   print "    Distance between $clusterJ and ancestral node = $branch2\n";

   for ($i = 0; $i < $numClusters; $i++){
        if (($arrayClusters[$i] ne $clusterI)
           && ($arrayClusters[$i] ne $clusterJ)){
           $d1 = $hashClusters{$arrayClusters[$i]}{$clusterI};
           $d2 = $hashClusters{$arrayClusters[$i]}{$clusterJ};
           # calculate new distance value
           $hashClusters{$merge}{$arrayClusters[$i]} =
                  ($d1 + $d2 - $hashClusters{$clusterI}{$clusterJ}) / 2;
           $hashClusters{$arrayClusters[$i]}{$merge} = $hashClusters{$merge}{$arrayClusters[$i]};
       }
   }

   $hashNewick{$merge} = "(".$hashNewick{$clusterI}.",".$hashNewick{$clusterJ}.")";
   # eliminate hash keys $clusterI and $clusterJ, decrement numClusters
   for ($i=0; $i<$numClusters;$i++){
     delete $hashClusters{$clusterI};
     delete $hashClusters{$clusterJ};
   }
   delete $hashNewick{$clusterI};
   delete $hashNewick{$clusterJ};
   
   # create new array of cluster values
   @arrayClusters = keys(%hashClusters);
   
   $numClusters --;
}

#print information
print "Distance between remaining clusters: ";
print $hashClusters{$arrayClusters[0]}{$arrayClusters[1]},
      $hashClusters{$arrayClusters[1]}{$arrayClusters[0]};
# should have two clusters remaining
print "\nNewick format: ";
for ($j = 0; $j < $numClusters; $j++){
   print $hashNewick{$arrayClusters[$j]}, " ";
}
