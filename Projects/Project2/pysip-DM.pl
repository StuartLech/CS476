#Exploring Bioinformatics: A Project-Based Approach
#Caroline St.Clair and Jonathan Visick

#Chapter 7 - Skills 1-5* 
# This implementation only implements UPGMA and only reads
# in distance matrices. 
# To make the Newick format easier to read, this implementation
# uses the first character in the fasta line of the sequences
# file as sequence names

# Modified significantly by Gunes 
# to allow distance matrices to be entered
# rather than distances to be computed!

%hashOriginal = (); # hash of hash to hold initial distance matrix value

%hashClusters = (); # hash of hash to hold matrix value
                    #    key = code for sequence 1
                    #    value = hash of sequence 2 code and distance

%hashNewick = (); # hash table to hold Newick code of cluster

%codehash = (); # this is the "inverse" of the codes, yielding back the original index

my($data); # temp variable to hold current data
my(@codes); # This is the character code for each original OTUs
my($numOTU); # This is the number of OTUs (or the length of the matrix)
my(@distmatrix); # This is the "matrix" version of the distance matrix!
my(@temp); # temporary array

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
    $codehash{$codes[$i]} = $i;
    $hashNewick{$codes[$i]} = $codes[$i];
    $hashClusters{$codes[$i]}{0} = $codes[$i];
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
	$hashOriginal{$codes[$i]}{$codes[$j]} = $distmatrix[$i][$j];
    }
    print "\n";
}

#########################################################
# Agglomerative Clustering Algorithm Implementing UPGMA
#     - use arithmetic distance if single items in clusters
#       otherwise uses UPGMA
#     - merge clusters with smallest distance
#     - repeat process until one cluster remains
#########################################################
$MAXDIST = 999999;
$numClusters = $numOTU;
while ($numClusters > 1){
   # Calculate distances and find smallest
   $smallest = $MAXDIST; $smallestI = 0; $smallestJ = 0;
   @arrayClusters = keys(%hashClusters);
   for ($i = 0; $i < $numClusters-1; $i++){
       for ($j = $i+1; $j < $numClusters; $j++){
	   # Debugging:
	   print "\nCluster $i: $arrayClusters[$i] Cluster $j: $arrayClusters[$j]\n";
	   $tempdist = 0;
	   # The following code determines distance between
	   # clusters as the Unweighted Average of all distances.
	   for ($k = 0; $k < length $arrayClusters[$i]; $k++){
	       for ($m = 0; $m < length $arrayClusters[$j]; $m++){
		   $tempdist = $tempdist + 
		       $hashOriginal{$hashClusters{$arrayClusters[$i]}{$k}}{$hashClusters{$arrayClusters[$j]}{$m}};
	       }
	   }    
	   $tempdist = $tempdist / ((length $arrayClusters[$i]) *
				    (length $arrayClusters[$j]));
	   if ($tempdist < $smallest){
	       $smallest = $tempdist;
	       $smallestI = $i;
	       $smallestJ = $j;
	   }
       }
   }
   # merge smallestI and smallestJ into clusters and add
   # to hash tables
   $clusterI = $arrayClusters[$smallestI];
   $clusterJ = $arrayClusters[$smallestJ];
   $merge = $clusterI.$clusterJ;
   print "Merging Clusters: $clusterI and $clusterJ with distance $smallest\n";
   $i=0;
   for ($j = 0; $j < length $clusterI; $j++){
       $hashClusters{$merge}{$i} = $hashClusters{$clusterI}{$j};
       $i++;
   }
   for ($j = 0; $j < length $clusterJ; $j++){
       $hashClusters{$merge}{$i} = $hashClusters{$clusterJ}{$j};
       $i++;
   }
   $hashNewick{$merge} = "(".$hashNewick{$clusterI}.",".$hashNewick{$clusterJ}.")";
   # eliminate old clusters, decrement numClusters
   delete $hashClusters{$clusterI};
   delete $hashClusters{$clusterJ};
   delete $hashNewick{$clusterI};
   delete $hashNewick{$clusterJ};
   $numClusters --;
}

#print hash structures
@arrayClusters = keys(%hashClusters);  #should only be one cluster left
#print "Cluster ", $arrayClusters[0], " values: ";
#for ($j = 0; $j < length $arrayClusters[0]; $j++){
#   print $hashClusters{$arrayClusters[0]}{$j}, " ";
#}
print "\nNewick format: ", $hashNewick{$arrayClusters[0]}, "\n";
