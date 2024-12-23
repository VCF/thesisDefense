#!/usr/local/bin/perl
BEGIN {
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser. From Steve's Primer3 program
    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
    push @INC, "/usr/local/etc/httpd/cgi-bin/pagelab/charles/Probability";
    push @INC, "/usr/local/etc/httpd/cgi-bin/pagelab/charles/CalcTheta";
}

# Staffa path for this file is:
# /usr/local/etc/httpd/cgi-bin/pagelab/charles/
# SimulateRH.cgi
# A Radiation Hybrid simulator with predictive capacity
# Charles Tilford, 2000
# Whitehead Institute / Massachusetts Institute of Technology

# My stupid reminders on how to use XS:
# Move to directory with .xs file and execute	: perl Makefile.PL
# Then execute									: make
# Find the #!@%@ .so file. It was burried in	: blib/arch/auto/Probability
# Move a copy to the Probability level, or put the above in @INC. That's it!
 
strict;
use CGI;
use GD;
use Probability;
use CalcTheta;
require "ParseQuery.pl";
# $cDebugger = 1;

&DEFINECONSTANTS;
&PARSE;


&FILEDUMP($Dump) if ($Dump);


if ($filename) {
	if ($noStore) {
		&LOADDATA;
		&DYNAMICLOAD;
	} else {
		&EXTERNALLOAD;
	}
}

&HELPME			if ($helpMe);
&STARTHTML;
&LOADDATA;
&FINDLINKGROUPS	if ($lodLink);
&SHOWSTATS		if ($ShowStats);
&SHOWVECTORS	if ($ShowVectors);
&PRINTOPTIONS;
&SEQCOMPARE		if ($seqCompare);
# $a = Probability::twoVector($myVectors{"sY1"}, $myVectors{"sY2"}); @t = split "\t", $a; foreach $i (@t) { print "$i<BR>\n"; }

if ($filename ne "") {
	print "<FONT COLOR=green>Your file was loaded. You may now use the options above to analyze it.</FONT><BR>\n";
	exit;
}
unless ($doPredict eq "" || $doPredict eq "No Prediction") {
	my $targ = \@vecOrder;
	@tempOrder = ();
	if ($myMarkers) {
		@tempOrder = split /[\r\n]/, $myMarkers;
		$targ = \@tempOrder;
	}
	@holding = (); my @unfound = ();
	for my $i (0..$#{$targ}) {
		next if (${$targ}[$i] eq "");
		if (exists $myVectors{ ${$targ}[$i] }) {
			push @holding, ${$targ}[$i];
		} else {
			push @unfound, ${$targ}[$i];
		}
	}
	if ($#unfound > -1) {
		my $tot = $#holding +1; my $unf = $#unfound + 1;
		print qq(<FONT COLOR=red><B>Could not find the following $unf markers ($tot were found) in the local vector file:</B>\n);
		print "<PRE>";
		for my $i (0..$#unfound) {
			print "\n" unless ($i % 6);
			printf ("%15.15s ", $unfound[$i]);
		}
		print "</FONT></PRE>";
	}
	$targ = \@holding;
	if ($doPredict =~ /all/i) {							# Do a full grid analysis
		if ($gridOrder =~ /alph/i) {					# Arrange markers alphabetically
			@tempOrder = sort {$a cmp $b} @{$targ};
			$targ = \@tempOrder;
		}
		if ($findLinkage) {
			&LINKAGETEST($targ, "PredictionOutput");
		} else {
			&GRIDPREDICT($targ, "PredictionOutput");
		}
	} elsif ($doPredict eq "Neighbors") {
		&NEIGHBORPREDICT($targ, "PredictionOutput");
	}
	exit;
}
if ($thePanel) {
	if ($scanStep > 0) {
		&STEPRUN($thePanel, "A","B");
	}
}


# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# This routine was made to verify that probabilities needed to be modified by N choose K options
# (they do - the sum of the probs is not one without them).
sub PROBTEST {
	my $hybs = 5;
	my $pa = 0.1;
	my $pb = 0.22;
	my $pc = 0.2;
	my $pd = 1 - $pa - $pb - $pc;
	my @fact;
	for $i (0..$hybs) {
		$fact[$i] = 1;
		for $j (1..$i) {
			$fact[$i] *= $j;
		}
	}
	for $a (0..$hybs) {
		my $Pa = Probability::simplePower($pa, $a);
		for $b (0..($hybs-$a)) {
			my $Pb = Probability::simplePower($pb, $b);
			for $c (0..$hybs-$a-$b) {
				my $Pc = Probability::simplePower($pc, $c);
				for $d (0..$hybs-$a-$b-$c) {
					next if ($a+$b+$c+$d != $hybs);
					my $Pd = Probability::simplePower($pd, $d);
					$P = $Pa * $Pb * $Pc * $Pd;
					$P *= $fact[$hybs] / ($fact[$a] * $fact[$b] * $fact[$c] * $fact[$d]);
					printf("%d%d%d%d = %f<BR>\n", $a, $b, $c, $d, $P);
					$TotP += $P;
				}
			}
		}
	}
	print "<BR>TOTAL = $TotP<BR>\n";
	die
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub FINDLINKGROUPS {
	my $minLOD		= shift;		# minimum LOD to form groupings
	my $markerRef	= shift;		# Reference to list of markers to be analyzed
	my $mode		= shift;		# Type of LOD calculation to use - "jones" specifies the standard method, otherwise the simulator is used
	
	$minLOD = $lodLink unless ($minLOD);
	$mode = $linkMethod unless ($mode);
	unless ($markerRef) {
		if ($myMarkers) {
			$myMarkers =~ s/ +/ /g;
			$markerRef = [split /[\r\n ]/, $myMarkers];		# Custom markers
		} else {
			$markerRef = \@vecOrder;
		}
	}
	
	my @markers = ();
	my @rejects = ();
	my %inThere = ();
	foreach $sts (@{$markerRef}) {
		next if (exists $inThere{ $sts });	# Prevent duplicates
		if (exists $myVectors{ $sts }) {	# Only consider RH typed markers
			push @markers, $sts;
			$inThere{$sts} = "Y";
		} else {
			push @rejects, $sts if (length($sts) > 0);
		}
		
	}
	
	my %group = ();					# Will track which group a marker is in
	my @freeGroups = ();			# The available free groups
	my $groupCounter = 1;			# Current open group
	my @theGroups = ();				# The actual groups of markers
	my @unlinked = ();				# Pool of totally unlinked markers
	my $totAnalyzed = 0;			# Number of pairwise comparisons
	my $tim = time;
	for $r (0..($#markers-1)) {
		my $row = $markers[$r];
		for $c (($r+1)..$#markers) {
			my $col = $markers[$c];
			my $Lod;
			if ($mode =~ /jones/i) {
				$Lod = CalcTheta::findLOD($myVectors{$row}, $myVectors{$col});
			} else {
				($Lod, $devLod) = split "\t", Probability::twoVector($myVectors{$row}, $myVectors{$col});
			}
			$Lod = int(0.5 + 1000 * $Lod) / 1000;
#			printf ("%10.10s vs. %10.10s = %s (%.3f)<BR>\n", $row, $col, $Lod, CalcTheta::findLOD($myVectors{$col}, $myVectors{$row}) );
			next if ($Lod < $minLOD);
			
			if (exists $group{$col} && exists $group{$row}) {	# both are already in a group
				next if ($group{$col} == $group{$row});			# Oops, its the same group...
				push @{$theGroups[$group{$row}]}, @{$theGroups[$group{$col}]};
				push @freeGroups, $group{$col};					# Indicate that the col group may now be re-populated
#				foreach $n (@freeGroups) {print "$n, ";} print "<BR>";
				my $slaggedGroup = $group{$col};
				foreach $sts (@{$theGroups[$group{$col}]}) {
					$group{$sts} = $group{$row};				# Reset the group numbers for the column group
				}
				@{$theGroups[$slaggedGroup]} = ();				# Clear the col group
			} elsif (exists $group{$col}) {						# Col in group, Row not
				push @{$theGroups[$group{$col}]}, $row;			# Add it to the group...
				$group{$row} = $group{$col};					# ... and set its indicator
			} elsif (exists $group{$row}) {						# Row in group, Col not
				push @{$theGroups[$group{$row}]}, $col;			# Add it to the group...
				$group{$col} = $group{$row};					# ... and set its indicator
			} else {											# Neither are in a group, start a new one
				my $groupNum = $groupCounter;
				unless ($#freeGroups < 0) {
					$groupNum = pop @freeGroups;				# Re-use a depopulated group
				} else {
					$groupCounter++;							# No free groups, make a new one
				}
				@{$theGroups[$groupNum]} = ($row, $col);
				$group{$row} = $group{$col} = $groupNum;
			}
		}
		push @unlinked, $row unless (exists $group{$row});		# No linkage found, put it in unlinked group
		$totAnalyzed += $#markers - $r;
	}
	
	$tim = time - $tim;
	my $rate = $tim > 0 ? int($totAnalyzed / $tim) : "Lots";
	printf ("%d pairwise comparisons in %d seconds = %s per second.<BR>\n", $totAnalyzed, $tim, $rate);
	for $i (0..$#theGroups) {
		next if ($#{$theGroups[$i]} < 0);
		print "<B><FONT COLOR=green>Group $i</FONT></B>";
		&TABLEPRINT(\@{$theGroups[$i]}, 10, 8);
	}
	print "<FONT COLOR=red><B>Unlinked Markers</B></FONT>";
	&TABLEPRINT(\@unlinked, 10, 8);
	if ($#rejects > -1) {
		print "<FONT COLOR=red><B>Rejected Markers</B></FONT>";
		&TABLEPRINT(\@rejects, 10, 8);
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub TABLEPRINT {
	my $valRef		= shift;	# ref to array of values to print
	my $w			= shift;	# width of each column
	my $c			= shift;	# number of columns
	print "<PRE>";
	for $j (0..$#{$valRef}) {
		print "\n" unless ($j % $c);
		printf ("%${w}.${w}s", ${$valRef}[$j]);
	}
	print "</PRE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADORDER {
	my $file	= shift;		# filename
	open (FILE, "$file") or die "Failure to open $file: $!\n";
	my @skip = ("", "Boundary", "Gap");
	my $currentChunk = "null";
	my %chunks = ();
	while (<FILE>) {
		chomp ($_);
		next unless ($_);
		next if ($_ eq "Boundary");
		if ($_ =~ /#/) {		# Designated comment or chunk zone
			s/[# ]//g;
			if (/chunk/i) {
				$currentChunk = $_;
			}
			next;
		}
		push @{$chunks{$currentChunk}}, $_;
	}
	close FILE;
	return \%chunks;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SEQCOMPARE {
	print "<CENTER><TABLE><TR><TD BGCOLOR=yellow WIDTH=300><CENTER>";
	print "<FONT COLOR=brick SIZE=+4><B>Please do not disturb</B></FONT></TABLE></CENTER><BR>\n";

	my $dir = "$HOME/Data/ChrYMarkers";		# The directory where the clones lurk
	my @fileList = split "\n", `ls $dir`;
	open (FILE, ">$HOME/Data/PlotData") or die "SEQCOMPARE: Can't write PlotData: $!\n";
	print FILE "Marker A\tMarker B\tLOD\tDistance (kb)\tMarker A Sites\tMarker B Sites\tTheta\tJones LOD\n";
	close FILE;
	
	undef(%primerSites); undef(%fullSites); 
	foreach $f (@fileList) {
		open (FILE, "$dir/$f") or die "SEQCOMPARE: Can't read $dir/$f: $!\n";
		my ($chunkName, $chunkSize) = split "\t", <FILE>;
		while (<FILE>) {
			chomp;
			my @temp = split "\t", $_;
			my $sts = shift @temp;
			$primerSites{$sts}	+= 	shift @temp;
			$fullSites{$sts}	+= $#temp + 1;
		}
		close FILE;
	}
	my $RHoA_order = &LOADORDER("$HOME/Data/order");
	
#	foreach $s (@{${$RHoA_order}{'Chunk01'}}) {print "$s<BR>";}
	foreach $f (@fileList) {
		$_ = $f; /.*(\d\d)/; my $num = $1;
		print "<FONT COLOR=green SIZE=+2><B>$f</B>";
		my %positions = (); my @group = (); my %primerSites = ();
		open (FILE, "$dir/$f") or die "SEQCOMPARE: Can't read $dir/$f: $!\n";
		($chunkName, $chunkSize) = split "\t", <FILE>;
		printf (" %s @ %.3f Mb</FONT><BR>\n", $chunkName, $chunkSize / 1000000);
		while (<FILE>) {
			chomp;
			my @temp = split "\t", $_;
			my $sts = shift @temp;
			$primerSites{$sts} =	shift @temp;
			@{$positions{$sts}} = @temp;
		}
		close FILE;
		
		
		my %tempRef = %{$RHoA_order};
		&PLOTSTS(\%positions, $chunkSize, \@{$tempRef{"Chunk$num"}});
		next if ($noPlace);
		
		foreach $sts (keys %positions) {
			push @group, $sts if (exists $myVectors{ $sts });
		}
		@group = sort {$a cmp $b} @group;
		open (FILE, ">>$HOME/Data/PlotData") or die "SEQCOMPARE: Can't append to PlotData: $!\n";
		print "<PRE><i><FONT COLOR=brown>Now checking: </FONT></i>";
		for $col (0..($#group - 1)) {
			print "\n" unless ($col % 8);
			$aMarker = $group[$col];
			printf ("%12s ", $aMarker);
			for $row (($col+1)..$#group) {
				$bMarker = $group[$row];
				my $dist = 999999999999;
				foreach $a (@{$positions{$aMarker}}) {
					foreach $b (@{$positions{$bMarker}}) {
						my $theDist = abs($a - $b);
						$dist = $theDist if ($theDist < $dist);		# Find the closest pair of markers
					}
				}
				next if ($dist == 0);								# Zero distance markers are a pain to deal with when plotting on log axes
				my @Output = ();
				my $totTwo = $myVectors{$aMarker} . $myVectors{$bMarker};
				my @tmpA = ($totTwo =~ /(2)/g);						# Count the combined number of twos
				$totTwo = $#tmpA + 1;
				next if ($totTwo > $maxAmbig);						# Too many Ambiguous
				my ($theLOD, $devLOD) = split "\t", Probability::twoVector($myVectors{$aMarker}, $myVectors{$bMarker});
				next if (abs($theLOD) == 50);
				my $jonesLOD = int(0.5 + 1000 * CalcTheta::findLOD($myVectors{$aMarker}, $myVectors{$bMarker}))/1000;
				next if (($theta = CalcTheta::findTheta($myVectors{$aMarker}, $myVectors{$bMarker})) == 5);
				printf FILE ("%s\t%s\t%s\t%.3f\t%s\t%s\t%.3f\t%s\n", $aMarker, $bMarker, int(0.5 + 1000*$theLOD)/1000,
					$dist/1000, $fullSites{$aMarker}, $fullSites{$bMarker}, $theta, $jonesLOD);
			}
		}
		close FILE;
		print "</PRE>";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PLOTSTS {
	my $posRef		= shift;	# [0] reference to hash containing base pair position location for each marker
	my $totalSize	= shift;	# [1] total size of genomic sequence, bp
	my $orderRef	= shift;	# [2] ref to array of RH ordered markers
	
	
#	foreach $s (@{$orderRef}) { print "$s<BR>"; }

	my $bar			= 4;	# height of rectangle representing DNA
	my $tick		= 6;	# height of ticks connecting markers to DNA bar
	my $tickgap		= 2;	# space between tick and marker name
	my $maxChars	= 8+8;	# Maximum number of characters to print in a name
	my $Border		= 20;	# empty space on edges of DNA
	my $kbPerPix	= 3;	# Number of kilobases represented by one pixel
	my $connect		= 50;	# distance over which connecting lines are drawn

	my $bpPerPix	= $kbPerPix * 500;
	$barBottom		= $bar + $tick + $tickgap + $wT * $maxChars;
	
	$height			= $barBottom + $connect + $wT * $maxChars;
	$width			= $Border*2 + $totalSize/$bpPerPix;
	$gph			= new GD::Image($width,$height);
	$white			= $gph->colorAllocate(255,255,255);
	$almostWhite	= $gph->colorAllocate(255,255,254);
	$black			= $gph->colorAllocate(0,0,0);
	$blue			= $gph->colorAllocate(0,0,255);
	$orange			= $gph->colorAllocate(255,128,0);

	$red			= $gph->colorAllocate(255,0,0);
	$yellow			= $gph->colorAllocate(255,255,0);
	$green			= $gph->colorAllocate(0,128,0);
	$purple			= $gph->colorAllocate(255,0,128);
	$ltblue			= $gph->colorAllocate(128,128,255);
	$gray			= $gph->colorAllocate(150,150,150);
	
	my $halfHeight	= int ($hT/2);
	my $tickY2		= $barBottom - $bar - 1;
	my $tickY1		= $tickY2 - $tick;
	my $textY		= $tickY1 - $tickgap;

	$gph->filledRectangle($Border,$barBottom-$bar,$width-$Border,$barBottom,$orange);
	
	my @toSort = ();
	my %markPos = %{$posRef};
	my $skipped = my $notRH = 0;
	foreach $sts (keys %markPos) {
		unless ($plotAll) {
			unless (exists $myVectors{ $sts }) {
				$skipped += $#{$markPos{$sts}} + 1; $notRH += $#{$markPos{$sts}} + 1;
				next;
			}
		}
		foreach $p (@{$markPos{$sts}}) {
			push @toSort, sprintf("%12s\t%s", $p, $sts);
		}
		$col{$sts} = $#{$markPos{$sts}} > 0 ? $blue : $black;
	}
	@toSort = sort {$a cmp $b} @toSort;
	
	if ($orderRef =~ /seq/) {			# Make the RH map agree totally with the sequence
		my @newRhOrder = ();
		foreach $e (@toSort) {
			my ($pos, $sts) = split "\t", $e;
			if (exists $myVectors{$sts}) {
				push @newRhOrder, $sts;
			}
		}
		$orderRef = \@newRhOrder;
	}

	push @toSort, "99999999999999\tDummy";
	my $lastTex = my $lastPos = 0;
	my $minSpace = $hT + 1;
	my $sepPad = 20;
	undef (%HoA_seqPlot);		# Hash holding image location of markers
	for $i (0..$#toSort) {
		my ($pos, $sts) = split "\t", $toSort[$i];
		my $pix = $Border + $pos / $bpPerPix;
		push @{$HoA_seqPlot{$sts}}, $pix;
		my $tex = $pix - $halfHeight;
		my $skew = 0;
		my ($posNext, $d) = split "\t", $toSort[$i+1];
		my $nextTex = $Border + ($posNext / $bpPerPix) - $halfHeight;
		if ($geneAlias{$sts} ne "") {
			$gph->stringUp(gdTinyFont,$tex,$textY - $wT*9,$geneAlias{$sts},$orange);
		}
		if ($lastTex + $minSpace > $tex) {		# Will the item conflict with the prevous marker?
			if ($lastTex + $minSpace * 2 < $nextTex) {
				$tex = $lastTex + $minSpace;	# We can squeeze the marker in between the flanking markers
			} else {
				$skipped++;
				next;							# Can't fit it in, skip it
			}
		} elsif ($tex + $minSpace > $nextTex) {	# There *will be* a conflict with the next marker
			if ($lastTex + $minSpace * 2 < $nextTex) {
				$tex = $nextTex - $minSpace;	# We can squeeze the marker in between the flanking markers
			}
		}
		$gph->line($tex+$halfHeight, $tickY1, $pix, $tickY2, $orange);
		$gph->stringUp(gdTinyFont,$tex,$textY,$sts,$col{$sts});
		my $gapDist = $pos - $lastPos;
		my $gapTag = int(0.5 + $gapDist / 1000) . " kb";
		my $tagTex = $Border + ($lastPos + $gapDist/2)/$bpPerPix - $wT * length($gapTag)/2;
		unless ($tagTex < $lastTex + $sepPad|| $gapDist  > 100000000) {			# Don't draw the gap size if it will be too close to the last marker 
			$gph->string(gdTinyFont,$tagTex,$textY,$gapTag,$red);
		}
		$lastTex = $tex;
		$lastPos = $pos;
	}

	my $cR = 5300;
	my $totalDist = 0;				# Total distance of the contig
	my $gapDist = -100 * $cR * log(1 - 0.6);		
	my $numContig = 0;				# Number (index) of RH contigs in the area
	my @A_orderDistance;
	my @A_orderName;
	@{$A_orderDistance[$numContig]} = (0);
	@{$A_orderName[$numContig]} = ($$orderRef[0]);
	for $i (1..$#{$orderRef}) {
		my $left = $$orderRef[$i-1];
		my $right = $$orderRef[$i];
		my $theta = CalcTheta::findTheta($myVectors{$left}, $myVectors{$right});
		my $dist = $gapDist;
		if ($theta <= 0.6) {
			$dist = -100 * $cR * log(1 - $theta);		
		} else {
			$numContig++;
			$totalDist = $dist = 0;
		}
		$totalDist += $dist;
	#	print "$right, $totalDist: ";
		push @{$A_orderDistance[$numContig]}, $totalDist;
		push @{$A_orderName[$numContig]}, $right;
	}
	
	for $c (0..$#A_orderDistance) {
		my $elem = $#{$A_orderName[$c]};			# Number of elements in this contig
		my $farLeft = $A_orderName[$c][0];		# Leftmost marker
		my $farRight = $A_orderName[$c][$elem];	# rightmost marker
		my $lp = 0;								# leftmost pixel
		my $rp = 30;							# rightmost pixel
		
		#Determine the left boundary:
		unless (exists $HoA_seqPlot{$farLeft}) {		# The left marker is not present in the sequence
			my $lc = 0;
			until (exists $HoA_seqPlot{$A_orderName[$c][$lc]} || $lc >= $elem)
				{ $lc++; }
			my $newLeft = $A_orderName[$c][$lc];
			$lp = $HoA_seqPlot{$newLeft}[0];
			print "Left marker $farLeft not present in sequence, using $newLeft instead.<BR>\n";
		} elsif ($#{$HoA_seqPlot{$farLeft}} == 0) {		# One position, set the left pixel to it
			$lp = $HoA_seqPlot{$farLeft}[0];
		} else {
			print "Multiple hits on left marker $farLeft.<BR>\n";
			$lp = $HoA_seqPlot{$farLeft}[0];
			if ($farLeft =~ /SP57/) {
				$lp = $HoA_seqPlot{'SP53-03'}[1];
				foreach $ff (@{$HoA_seqPlot{$farLeft}}) {
					print "$farLeft at $ff<BR>";
				}
			}
		}

		
		unless (exists $HoA_seqPlot{$farRight}) {		# The right marker is not present in the sequence
			my $rc = $elem;
			until (exists $HoA_seqPlot{$A_orderName[$c][$rc]} || $rc <= 0)
				{ $rc--; }
			my $newRight = $A_orderName[$c][$rc];
			$rp = $HoA_seqPlot{$newRight}[0];
			print "Right marker $farRight not present in sequence, using $newRight instead.<BR>\n";
		} elsif ($#{$HoA_seqPlot{$farRight}} == 0) {	# One position, set the right pixel to it
			$rp = $HoA_seqPlot{$farRight}[0];
		} else {
			print "Multiple hits on right marker $farRight.<BR>\n";
			my @t = @{$HoA_seqPlot{$farRight}};
			$rp = $t[$#t];
			$rp = $t[0] if ($farRight =~ /sY707/);
			if ($farRight =~ /sY796/) {
				# $rp = $t[0];
				foreach $ff (@t) {
					print "$farLeft at $lp, $farRight at $ff<BR>";
				}
			}
		}
		
		# These are left boundaries that are being manually defined:
		my %leftBastards = (	'sY13' => 'SP05-23,0',		'sY62' => 'sY733,0',
								'SP22-18' => 'SP22-33B,0',	'SP39-28' => 'sY101,0',
								'SP19-07' => 'sY605,0',		'sY71' => 'sY605,0',
								'SP25-39' => 'sY82,0',		'xxx' => 'xxx,0',
								
								'sY554' => 'sY554,0',		'sY820' => 'sY820,0',
								'SP54-47' => 'SP54-47,0',	'sY166' => 'sY166,0',
								'sY698' => 'sY698,0',		'sY908' => 'sY908,0',		'sY844' => 'sY804,0',
								'SP57-46' => 'sY1077,1',	'sY903' => 'sY903,2',
								'xxx' => 'xxx,0',		'xxx' => 'xxx,0',
								'xxx' => 'xxx,0',		'xxx' => 'xxx,0',
								);
								
		foreach $lb (keys %leftBastards) {
			if ($farLeft =~ /$lb/) {
				my ($sts, $num) = split ',', $leftBastards{$lb};
				$lp = $HoA_seqPlot{$sts}[$num];
			}
		}
		# These are right boundaries that are being manually defined:
		my %rightBastards = (	'sY1029' => 'sY1011,0',		'sY715' => 'sY715,0',
								'sY71' => 'SP19-19,0',		'SP19-07' => 'SP19-19,0',	
								'SP37-38' => 'sY100,0',		'SP41-31' => 'SP41-34,0',
								'sY911' => 'SP36-04,0',		'SP48-21' => 'SP47-28B,0',
								'sY924' => 'sY941,0',		'sY927' => 'sY941,0',
								
								'sY820' => 'sY820,0',		'sY554' => 'sY579,0',
								'sY166' => 'sY166,0',		'SP54-47' => 'SP54-11,0',
								'SP57-36' => 'SP57-36,0',	'sY707' => 'sY904,1',		'sY844' => 'sY844,0',	'sY908' => 'sY908,0',
								'sY796' => 'sY796,1',
								'xxx' => 'xxx,0',		'xxx' => 'xxx,0',
								'xxx' => 'xxx,0',		'xxx' => 'xxx,0',
								);
								
		foreach $rb (keys %rightBastards) {
			if ($farRight =~ /$rb/) {
				my ($sts, $num) = split ',', $rightBastards{$rb};
				$rp = $HoA_seqPlot{$sts}[$num];
			}
		}

		my $seqPix = $rp - $lp;							# Number of pixels occupied by the sequence
		my $rhSize = $A_orderDistance[$c][$elem];
		my $f_scaler = 1;
		if ($rhSize != 0) {
			$f_scaler = $seqPix / $rhSize;
		}
		my $lowerTick = $barBottom + $connect;
		my $lowerText = $lowerTick + $tick + $tickGap;
		my $oldPix = 0;
		for $i (0..$elem) {
			my $pos = $A_orderDistance[$c][$i];
			my $pix = $lp + $f_scaler * $pos;			# The x-coord on the gif image
			my $tex = $pix - $halfHeight;
			
			$gph->line($pix, $lowerTick, $pix, $lowerTick+$tick, $orange);
			my $sts = $A_orderName[$c][$i];
		#	print "$sts ";
			my $len = $wT * length($sts);
			$gph->stringUp(gdTinyFont,$tex,$lowerText+$len,$sts,$black);
			if ($i > 0 && $pos - $A_orderDistance[$c][$i-1] < $gapDist) {
				$gph->filledRectangle($oldPix,$lowerTick-$bar,$pix,$lowerTick,$orange);
			}
			my $minDist = 999999999;
			my $seqPix = 0;
			foreach $sP (@{$HoA_seqPlot{$sts}}) {		# For multiple locations, find the one closest to the plotted value
				my $d = abs($pix - $sP);
	#			print "$d = abs($pix - $sP)<BR>" if ($sts =~ /200/);
				if ($d < $minDist) {
					$minDist = $d;
					$seqPix = $sP;
				}
			}
			$gph->line($seqPix, $barBottom, $pix, $lowerTick-$bar, $orange) if ($seqPix);
			$oldPix = $pix;
		}
	}

	
	$plotCounter++;
	my $fileName = "Sequence-$plotCounter-$$.gif";
	open (FILE, ">$DataLoc/$fileName") or die "Failure to write to $fileName: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$DataLoc/$fileName"><BR>);
	print "<i>$skipped markers were not drawn, $notRH because they have not been RH mapped.</i><BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - ref to the chromosome
# [1] - name of "MarkerA"
# [2] - name of "MarkerB"
sub STEPRUN {
	my $chromosome	= shift;
	my $MarkerA		= shift;
	my $MarkerB		= shift;
	my $fileN		= "SimulationOutput";
	if ($scanLow < 0 || $scanHigh < $scanLow or $scanStep <= 0) {
		print "<FONT COLOR=red>Uh, you seem to have put in flaky values for the scan settings. Please re-enter such that:<BR>\n";
		print "The low value is at least zero<BR>\n";
		print "The high value is greater than the low one<BR>\n";
		print "The step value is greater than zero.<BR></FONT>\n";
		return;
	}
	my @Dist = (); my %Values = ();
	print "<B>Distances being used:</B><BR>";
	if ($logDist) {
		my $minY = 1 - 1/exp($lambda * $scanLow);
		my $maxY = 1 - 1/exp($lambda * $scanHigh);
		my $scanSpace = ($maxY - $minY) / ($scanStep - 1);
		my $x = 0;
		for (my $y = $minY; $#Dist < $scanStep - 1; $y += $scanSpace) {		# Use $#Dist as primary counting mechanism
			$x =  abs(int (0.5 -log(1 - $y) / $lambda));					# abs helps avoid negative zero
			push @Dist, $x;
			print "$x kb<BR>";
		}
	} else {
		my $scanSpace = ($scanHigh - $scanLow) / ($scanStep - 1);
		for (my $i = $scanLow; $i <= $scanHigh; $i += $scanSpace) {
			push @Dist, $i;
			print "$i kb<BR>";
		}
	}
	my @checking = ();
	my $plotWidth	= 250;
	my $plotHeight	= 200;
	if ($chromosome eq "All") {
		foreach my $n (@theChr) {
#			print "--> $n<BR>";
			push @checking, $Chromosomes{$n};
		}
	} else {
		@checking = ($Chromosomes{$chromosome});
		$plotWidth	= 600;
		$plotHeight	= 400;
	}
	&SAVEPREP($fileN);
	$theTime = time;
	print "<TABLE BORDER=1>";	$cols = 3;
	for my $chr (0..$#checking) {
		next unless (defined $chr);
		print "<TR>" unless ($chr % $cols);
		print "\n<TD>";
		$checking[$chr]->calcProbs(\%GV, \@Dist, \%Values);
		my $ABequal = 1;
		foreach my $d (@Dist) {
			if ($Values{"A$nullToken"}[$d] != $Values{"B$nullToken"}[$d]) {
				$ABequal = 0;
				last;
			}
		}
		if ($ABequal) {
			delete $Values{"B$nullToken"};
			@{$Values{"A$nullToken & B$nullToken"}} = @{$Values{"A$nullToken"}};
			delete $Values{"A$nullToken"};
		}
		&SAVEDATA(\%Values, \@Dist, $checking[$chr], $fileN);
		my @plots = ( \%Values ); my @labels = ($checking[$chr]->{Name});
		&XYPLOT(\@labels, \@plots, \@Dist, $plotWidth, $plotHeight, "$DataLoc/Plot$chr$$.gif");
	}
	print "</TABLE>\n";
	$theTime = time - $theTime;
	$theOperations = ($#checking + 1) * ($#Dist+1);
	&SAVEEND($fileN);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SHOWSTATS {
	print "<TABLE BORDER=1>";
	my $cols = 4; my $count = 0; my %sameRF = ();
	for my $sts (0..$#Statistics) {
		print "<TR>" unless ($count % $cols);
		print "\n<TD>";
		my @plots = ( \%{$Statistics[$sts]->{Probability}} );
		$_ = $Statistics[$sts]->{Model}; s/\*/\-/g; s/\|/ /g;
		my $name = $_;
		my @labels = ( $name ); my @content = ();
		for (my $i = 0; $i <= length($name)-1; $i += 2) {	# Take every other character (should be marker name) and collect them
			push @content, substr($name, $i,1);
		}
		@content = sort {$a cmp $b} @content;
		$sig = join '', @content;
#		print "$name->$sig.";
		push @{$sameRF{$sig}}, $sts;
		&XYPLOT(\@labels, \@plots, \@{$Statistics[$sts]->{Distance}}, 200, 150, "$DataLoc/Plot$sts$$.gif");
		$count++;
	}
	foreach $sig (keys %sameRF) {
		next if ($#{$sameRF{$sig}} < 1);
		print "<TR>" unless ($count % $cols);
		print "\n<TD>";
		my @labels = ();
		my @plots = ();
		foreach $sts (@{$sameRF{$sig}}) {
			push @plots,	\%{$Statistics[$sts]->{Probability}};
			$_ = $Statistics[$sts]->{Model}; s/\*/\-/g; s/\|/ /g;
			push @labels,	$_;
#			print "$_, ";
		}
		&XYPLOT(\@labels, \@plots, \@{$Statistics[$sts]->{Distance}}, 200, 150, "$DataLoc/MultiPlot$sts$$.gif");
		$count++;
	}
	print "</TABLE>\n";
	exit;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SHOWVECTORS {
	my @vecKeys = @vecOrder;
	@vecKeys = @vecAlpha if ($AlphaSort);
	print "<H2><CENTER>Vectors stored locally:</CENTER></H2><PRE>";
	foreach $n (@vecKeys) {
		printf ("%15.15s ",$n);
		$_ = $myVectors{$n};
		my $swapOut = "<FONT COLOR=yellow>0</FONT>";
		s/0/$swapOut/g;
		my $swapOut = "<FONT COLOR=red>2</FONT>";
		s/2/$swapOut/g;
		print "$_\n";
	}
	print "</PRE>\n";
	if ($AlphaSort) {
		print qq(<FONT COLOR=orange><CENTER>View the list in the <A HREF="SimulateRH.cgi?ShowVectors=Y&dum=$$">original order</A> as found in file.</CENTER></FONT>\n);
	} else {
		print qq(<FONT COLOR=orange><CENTER>Sort the above list <A HREF="SimulateRH.cgi?ShowVectors=Y&AlphaSort=Y&dum=$$">alphabetically</A></CENTER></FONT>\n);
	}
	exit;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - filename
sub FILEDUMP {
	print "<html><TITLE>Text Dump for file $_[0]</TITLE><BODY BGCOLOR=White>\n";
	print "<PRE>";
	open (FILE, "$DataLoc/$_[0]") or die "Failure to open $DataLoc/$_[0]: $!\n";
	while (<FILE>) { print;}
	print "</PRE>";
	exit;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub NCHOOSEK {
	my $n	= shift;		# [0] - n, the total number of possibilities
	my $k	= shift;		# [1] - k, the number that you want to choose out of n
	my $res = 1;
	for my $i (1..$n) {
		$res *= $i;
	}
	for my $i (1..$k) {
		$res /= $i;
	}
	for my $i (1..($n-$k)) {
		$res /= $i;
	}
	return $res;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - filename where data should be saved
sub PREDICTPREP {
	my $saveFile	= shift;
	open (FILE, ">$DataLoc/$saveFile") or die "Failure to initialize $saveFile: $!\n";
	print FILE "START\n\n";
	print FILE "Value\tLOD Cutoff Limit\t$lodCutoff\n";
	print FILE "Value\tMaximum Ambiguous Results\t$maxAmbig\n";
	print FILE "Value\tMaximum Linkage Distance\t$maxLinkage\n";
	print FILE "Value\tNumber Next-Best Matches\t$keepNextBest\n";
	print FILE "\n";
	close FILE;	
	undef (%HitCount);
	undef (%lodSum1);
	undef (%lodSum2);
	undef (%distSum1);
	undef (%distSum2);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREDICTADD {
	my $saveFile	= shift;	# [0] - filename where data should be saved
	my $markerA		= shift;	# [1] - name of marker A 
	my $markerB		= shift;	# [2] - name of marker B 
	my $dataRef		= shift;	# [3] - ref to array holding data
	open (FILE, ">>$DataLoc/$saveFile") or die "PREDICTADD: Failure to append to $saveFile: $!\n";
	print FILE "$markerA\t$markerB";
	foreach $n (@{$dataRef}) {
		print FILE "\t$n";
	}
	print FILE "\n";
	close FILE;	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREPSTATVALS {
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# This routine prepares the statistics for random cost ordering.
# Because the matrix is non-symetrical, the average value ((X,Y) + (Y,X))/2 is taken for both X,Y and Y,X
# Only points where both X and Y have data are used - the remainder are set to a large value ($bogusHit) to avoid their contribution to ordering

# NEED TO MAKE FILES IN TMP FOLDERS, RATHER THAN IN MAIN ONE

sub STATOUT {
	my $saveFile	= shift;	# [0] - filename where data should be saved
	my $theOrder = \@allModels;
	my $statTime = time;
	$bogusHit = 45;
	
	# Calculate values for the purpose of random cost ordering
	for $a (0..($#allModels-1)) {
		my $topModel = $allModels[$a];
		$val{$topModel}{$topModel} = 0;
		for $b (($a+1)..$#allModels) {
			my $secModel = $allModels[$b];
			if (exists $HitCount{$topModel}{$secModel} && exists $HitCount{$secModel}{$topModel}) {
				my $hits	= $HitCount{$topModel}{$secModel} + $HitCount{$secModel}{$topModel};
				my $avgLod	= ($lodSum1{$topModel}{$secModel}  + $lodSum1{$secModel}{$topModel}) / $hits;
				$val{$topModel}{$secModel} = $val{$secModel}{$topModel} = sprintf("%.2f", $avgLod);
			} else {	# At least one of the models is screwy.
				$val{$topModel}{$secModel} = $val{$secModel}{$topModel} = $bogusHit;
			}
		}
	}
	my @havesData = ();		# Models that actually have useful information in them
	my @havenotsData = ();	# Models that only contain "$bogusHit" for values
	
	foreach $topModel (@{$theOrder}) {
		my $hit = 0;
		foreach $secModel (@{$theOrder}) {
			next if ($topModel eq $secModel);
			if ($val{$topModel}{$secModel} != $bogusHit) {
				$hit = 1;
				last;
			}
		}
		if ($hit != 0) {
			push @havesData, $topModel;
		} else {
			push @havenotsData, $topModel;
		}
	}
	
	open (FILE, ">$HOME/RCstats") or die "STATOUT: Failure to write to $HOME/RCstats: $!\n";
	foreach $topModel (@havesData) {
		print FILE "$topModel\t";
		foreach $secModel (@havesData) {
			print FILE "$val{$topModel}{$secModel}\t";
		}
		print FILE "\n";
	}
	close FILE;
	
	print "<FONT COLOR=green>Data written to file RCstats\n";
	my $t = time;
	my $dummyOut = (`/usr/local/etc/httpd/cgi-bin/pagelab/wu_data/src/rand_cost < $HOME/RCstats > $HOME/RCneighbors`);
	$t = time - $t;
	print ", and models ordered by random cost algorithm in $t seconds.</FONT><BR>\n";
	
	&LOADNEIGHBORS;

	for my $y (0..$#neighbors) {
		my $topModel = $neighbors[$y];
		for my $x (0..$#neighbors) {
			my $secModel = $neighbors[$x];
			if ($HitCount{$topModel}{$secModel} != 0) {
				my $hits	= $HitCount{$topModel}{$secModel};
				my $avgLod	= $lodSum1{$topModel}{$secModel} / $hits;
				my $avgDst	= $distSum1{$topModel}{$secModel} / $hits;
				my $devIn	= $avgLod*$avgLod - 2 * $avgLod *  $lodSum1{$topModel}{$secModel} / $hits + $lodSum2{$topModel}{$secModel} / $hits;
				my $devLod	= sqrt(int(0.5 + 1000000*$devIn)/1000000);
				my $devIn	= $avgDst*$avgDst - 2 * $avgDst * $distSum1{$topModel}{$secModel} / $hits + $distSum2{$topModel}{$secModel} / $hits;
				my $devDist	= sqrt(int(0.5 + 1000000*$devIn)/1000000);
				$AverageLODs{$topModel}{$secModel}		= $avgLod;
				$AverageDist{$topModel}{$secModel}		= $avgDst;
				$AverageLODsSD{$topModel}{$secModel}	= $devLod;
				$AverageDistSD{$topModel}{$secModel}	= $devDist;
			} else {
				$AverageLODs{$topModel}{$secModel}		= "";
				$AverageDist{$topModel}{$secModel}		= "";
				$AverageLODsSD{$topModel}{$secModel}	= "";
				$AverageDistSD{$topModel}{$secModel}	= "";
			}
		}
	}
	print "<FONT COLOR=green>Average values and standard deviations calculated for LODs and distances\n";
	my @Variables = ("AverageLODs", "AverageDist", "AverageLODsSD", "AverageDistSD");
	open (FILE, ">$HOME/StatisticTables") or die "STATOUT: Failure to write to $HOME/StatisticTables: $!\n";
	foreach my $v (@Variables) {
		print FILE "Variable\t$v\nOrder";
		for my $x (0..$#neighbors) {
			print FILE "\t$neighbors[$x]";
		}
		print FILE "\n";
		for my $y (0..$#neighbors) {
			my $topModel = $neighbors[$y];
			print FILE "$topModel";
			for my $x (0..$#neighbors) {
				my $secModel = $neighbors[$x];
				my $val = ${$v}{$topModel}{$secModel};
				if ($val ne "") {
					printf FILE ("\t%.2f", ${$v}{$topModel}{$secModel});
				} else {
					print FILE "\t";
				}
			}
			print FILE "\n";
		}
		print FILE "\n";
	}
	close FILE;
	print ", and written to file StatisticTables.</FONT><BR>\n";
	&DRAWAVGLOD;
	$statTime = time - $statTime;
	print "<FONT COLOR=brown>Total time devoted to statistics calculation and display was $statTime seconds.</FONT>\n";
	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADSTATS {
	open (FILE, "$HOME/StatisticTables") or die "STATOUT: Failure to load $HOME/StatisticTables: $!\n";
	my $theVar = "";			# The current variable being considered.
	my @theOrder = ();			# The order of models for the variable
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		my @t = split "\t", $_;
		if ($t[0] eq "Variable") {
			$theVar = $t[1];
			%{$theVar} = ();
			next;
		}
		if ($t[0] eq "Order") {
			$dum = shift @t;
			@theOrder = @t;
			next;
		}
		my $topModel = shift @t;
		for $y (0..$#theOrder) {
			${$theVar}{$topModel}{$theOrder[$y]} = $t[$y];
		}
	}
	close FILE;
	print "<FONT COLOR=green>All model statistics succefully loaded.<BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADNEIGHBORS {
	@neighbors = @origNeighbors = (); my %NeighAdd = ();
	open (FILE, "$HOME/RCneighbors") or die "STATOUT: Failure to read $HOME/RCneighbors: $!\n";
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		my $mod = $_;
		my $l = length($mod);
		$mod = substr($mod,0,$l-1) if (($l+1) % 2 != 0);	# The first model is being trashed for some reason... an extra odd character is being added to the name. I think rand_cost is doing this...
		push @neighbors, $mod;
		push @origNeighbors, $mod;
		$NeighAdd{$mod} = 1;
	}
	foreach $n (@allModels) {			# Add the unsorted neighbors to the end
		push @neighbors, $n if ($NeighAdd{$n} != 1);
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DRAWAVGLOD {
	my $saveFile	= shift;	# [0] - filename where data should be saved
	&LOADNEIGHBORS;
	&LOADSTATS;
	my $xOrder = \@neighbors;
	my $yOrder = \@origNeighbors;
	&STATPLOT($yOrder, $xOrder);
	for $y (0..$#{$yOrder}) {
		$topModel = ${$yOrder}[$y];
		my $yl = $y * $hT + $Border;
		my $yr = $yl + $hT - 1;
		for $x (0..$#{$xOrder}) {
			$secModel = ${$xOrder}[$x];
			my $xl = $x * $hT + $Border;
			my $xr = $xl + $hT - 1;
			my $avgLod = $AverageLODs{$topModel}{$secModel};
			if ($y == $x) {
				$gph->filledRectangle($xl,$yl,$xr,$yr,$white);
			} elsif ($avgLod ne "") {
				$gph->filledRectangle($xl,$yl,$xr,$yr,$LODcols[ int(0.5 + $avgLod / $stepFactor) ]);
			} else {		# No data for this point
				$gph->filledRectangle($xl,$yl,$xr,$yr,$ltblue);
			}
		}
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub STATPLOT {
	my $yMarkers	= shift;
	my $xMarkers	= shift;
	$edge			= 30;
	$Border			= $wT * $edge;
	$height			= $hT * ($#{$yMarkers}+1) + $Border*2;
	$width			= $hT * ($#{$xMarkers}+1) + $Border*2;
	$gph			= new GD::Image($width+$fudge,$height);
	$white			= $gph->colorAllocate(255,255,255);
	$black			= $gph->colorAllocate(0,0,0);
	$red			= $gph->colorAllocate(255,0,0);
	$orange			= $gph->colorAllocate(255,128,0);
	$yellow			= $gph->colorAllocate(255,255,0);
	$green			= $gph->colorAllocate(0,128,0);
	$blue			= $gph->colorAllocate(0,0,255);
	$purple			= $gph->colorAllocate(255,0,128);
	$ltblue			= $gph->colorAllocate(128,128,255);
	$gray			= $gph->colorAllocate(150,150,150);
	
	for my $i (0..$#{$yMarkers}) {
		my $n			= ${$yMarkers}[$i];
		my $l			= length($n);
		my $polarity	= (substr($n,0,1) eq "B") ? 1 : 0;		# Decide if model begins with "A" (0) or "B" (1)
		my $y 			= $i * $hT + $Border;
		my $x 			= $Border - 2 - $wT * length($n);
		for $letter (0..($l-1)) {
			my $delta = $letter * $wT;
			my $col = $black;
			my $c = substr($n, $letter, 1);
			$c = substr($n, $l-1-$letter, 1) if ($polarity);	# Read the model backwards
			if ($c eq 'B') {
				$col = $blue;
			} elsif ($c eq '|') {
				$col = $red;
			} elsif ($c eq '*') {
				$col = $green;
			}
			$gph->string(gdTinyFont,$x+$delta				,$y,$c,$col);
			$gph->string(gdTinyFont,$width-$Border+2+$delta	,$y,$c,$col);
		}
	}
	for my $i (0..$#{$xMarkers}) {
		my $n			= ${$xMarkers}[$i];
		my $l			= length($n);
		my $polarity	= (substr($n,0,1) eq "B") ? 1 : 0;		# Decide if model begins with "A" (0) or "B" (1)
		my $x			= $i * $hT + $Border;
		my $y			= $height - $Border + 2 + $wT * $l;
		for $letter (0..($l-1)) {
			my $delta = $letter * $wT;
			my $col = $black;
			my $c = substr($n, $letter, 1);
			$c = substr($n, $l-1-$letter, 1) if ($polarity);	# Read the model backwards
			if ($c eq 'B') {
				$col = $blue;
			} elsif ($c eq '|') {
				$col = $red;
			} elsif ($c eq '*') {
				$col = $green;
			}
			$gph->stringUp(gdTinyFont,$x,$Border-2 - $delta	,$c,$col);
			$gph->stringUp(gdTinyFont,$x,$y - $delta		,$c,$col);
		}
	}
	@LODcols = ();
	$maxPoor = 3; $stepFactor = 0.1;
	for ($p = 0; $p <= ($maxPoor+.01); $p += $stepFactor) {
		my $index = int(0.5  + $p / $stepFactor);
		$greenComponent = int (0.5 + 255 * $p / $maxPoor);			# Will make a gradient from red to yellow
		$newCol[$index] = $gph->colorAllocate(255,$greenComponent,0);
		$LODcols[$index] = $newCol[$index];
#		print "$index ";
	}
	$minGray = 3; $maxGray = 50;
	for $p ($minGray..$maxGray) {
		my $intensity = int (128 - 128 * ($p - $minGray) / ($maxGray - $minGray));
#		print "$p = $intensity<BR>";
		$grayCol[$p] = $gph->colorAllocate($intensity,$intensity,$intensity);
	}
	for ($p = $maxPoor + $stepFactor; $p <= 100; $p += $stepFactor) {
		my $index = int(0.5  + $p / $stepFactor);
		$gind = ($p > $maxGray) ? $maxGray : int(0.5 + $p);
#		print "$index ";
		$LODcols[$index] = $grayCol[$gind];
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREDICTSTAT {
	my $dataRef		= shift;	# [0] - ref to array holding data
	my ($lodBest, $lB, $distBest, $dB, $modelBest) = split "!!", ${$dataRef}[0];
	for my $n (1..$#{$dataRef}) {
		my ($lod, $l, $dist, $d, $model) = split "!!", ${$dataRef}[$n];
		$HitCount{$modelBest}{$model}++;
		my $lodDelta = $lod  - $lodBest;
		$lodSum1{$modelBest}{$model}	+= $lodDelta;
		$lodSum2{$modelBest}{$model}	+= $lodDelta * $lodDelta;
		my $distDelta = $dist - $distBest;
		$distSum1{$modelBest}{$model}	+= $distDelta;
		$distSum2{$modelBest}{$model}	+= $distDelta * $distDelta;
	}	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREDICTEND {
	my $saveFile	= shift;	# [0] - filename where data should be saved
	my $otherStuff	= shift;	# [1] - additional text to append after the STOP argument
	open (FILE, ">>$DataLoc/$saveFile") or die "PREDICTEND: Failure to append to $saveFile: $!\n";
	print FILE "\n\nSTOP\n";
	print FILE "\n$otherStuff\n";
	close FILE;	
	my $theSize = -s "$DataLoc/$saveFile";
	$theSize = int (0.5 + 100 * $theSize/(1024*1024))/ 100;
	print qq(<FONT COLOR=blue>Data saved to file <A HREF="SimulateRH.cgi?Dump=$saveFile&dum=$$" target=_blank>$saveFile</A>, total size $theSize Mb.</FONT>\n);
	print " <FONT COLOR=red SIZE=-1>This is far too big to view. Drag link to desktop to copy it to your hard disk.</FONT>\n" if ($theSize > 10);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREDICTINIT {
	my $numComp		= shift;	# [0] - number of comparisons that will be done
	my $saveFile	= shift;	# [1] - filename where data should be saved
	$pTime 			= $lastTime = time;
	$totalVariants = 0; $uniqueVariants = 0;
	my $units;
	my $runTime = $numComp * $Preferences{avgMarkerSpeed};
	if ($runTime > 15) {
		my $fin = time + $runTime;
		print "<CENTER><FONT COLOR=blue SIZE=+1><B>Simulation run time estimated at ";
		&CALCTIME(\$runTime, \$units);
		print "$runTime $units.<BR>\n";
		if ($units !~ /sec/i) {
			print "<TABLE><TR><TD BGCOLOR=yellow WIDTH=300><CENTER>";
			print "<FONT COLOR=brick SIZE=+4><B>Please do not disturb</B></FONT></TABLE><BR>\n";
			my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($fin);
			my @days = ('Sunday', 'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday');
			my $phase = "am";
			if ($hour > 12) { $phase = "pm"; $hour -= 12; }
			$min = "0" . $min if ($min < 10);
			print "Estimated finish at ${hour}:${min} $phase on $days[$wday]<BR>\n";
		}
		print "</FONT></CENTER></B>\n";
	}
	&PREDICTPREP($saveFile);
}

# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub GRIDPLOT {
	my $toTest		= shift;	# [0] - ref to array of known vectors that will be tested
	$edge			= 12;
	$Border			= $wT * $edge;
	$height			= $hT * ($#{$toTest}+1) + $Border*2;
	$width			= $hT * ($#{$toTest}+1) + $Border*2;
	$gph			= new GD::Image($width+$fudge,$height);
	$white			= $gph->colorAllocate(255,255,255);
	$black			= $gph->colorAllocate(0,0,0);
	$red			= $gph->colorAllocate(255,0,0);
	$green			= $gph->colorAllocate(0,128,0);
	$blue			= $gph->colorAllocate(0,0,255);
	$orange			= $gph->colorAllocate(255,128,0);
	$purple			= $gph->colorAllocate(255,0,128);
	$gray			= $gph->colorAllocate(150,150,150);
	
	open (FILE,"$HOME/Spectrum.gif") || die qq(Can't find "Spectrum.gif": $!\n);
	my $spect = newFromGif GD::Image(FILE) || die "Can't generate GIF from Spectrum.gif: $!\n";
	close FILE;
	@col = ();
	
	$numDivisions = 50;
	$colDiv = int ($maxLinkage / $numDivisions);
	for $i (0..($numDivisions+1)) {
		$index = $spect->getPixel(int(200 * $i/$numDivisions),1);
		$col[$i] = $gph->colorAllocate($spect->rgb($index));
	}
	my $step = 50;
	for (my $i=0; $i <= $maxLinkage / $step; $i++) {
		my $y = $i * $hT;
		my $n = $i * $step;
		my $c = $n / $colDiv;
		$n = " " . $n . " kb";
		$gph->filledRectangle($0,$y,$hT-2,$y+$hT-2,$col[ $c ]);
		$gph->string(gdTinyFont,$hT,$y,$n,$black);
	}
	push @col, $gray; 	
	%markPos		= ();
	for my $i (0..$#{$toTest}) {
		my $n			= ${$toTest}[$i];
		$markPos{ $n }	= $i;
		my $x			= $i * $hT + $Border;
		my $y			= $height - $Border + 2 + $wT * length($n);
		$gph->stringUp(gdTinyFont,$x,$Border-2	,$n,$black);
		$gph->stringUp(gdTinyFont,$x,$y			,$n,$black);
		$y 				= $x;
		$x 				= $Border - 2 - $wT * length($n);
		$gph->string(gdTinyFont,$x				,$y,$n,$black);
		$gph->string(gdTinyFont,$width-$Border+2,$y,$n,$black);
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LINKAGETEST {
	my $toTest		= shift;	# [0] - ref to array of known vectors that will be tested
	my $saveFile	= shift;	# [1] - filename where data should be saved
	my $totMark		= $#{$toTest} +1;
	my $numMrk = ($totMark - 1) * $totMark / 2;
	&PREDICTINIT($numMrk, $saveFile);

	my $edge		= 12;
	my $Border		= $wT * $edge;
	my $height		= $hT * ($#{$toTest}+1) + $Border*2;
	my $width		= $hT * ($#{$toTest}+1) + $Border*2;
	my $gph			= new GD::Image($width+$fudge,$height);
	my $white		= $gph->colorAllocate(255,255,255);
	my $black		= $gph->colorAllocate(0,0,0);
	my $red			= $gph->colorAllocate(255,0,0);
	my $pink		= $gph->colorAllocate(255,128,128);
	my $green		= $gph->colorAllocate(0,128,0);
	my $blue		= $gph->colorAllocate(0,0,255);
	my $purple		= $gph->colorAllocate(255,0,128);
	my $gray		= $gph->colorAllocate(150,150,150);

	my @linkCols = my @unlinkCols = (); %drawCols = ();
	my $maxPoor = 3; my $stepFactor = 0.1;
	if ($fancyColor) {
		for ($p = 0; $p <= ($maxPoor+.01); $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$in = int (0.5 + 128 * $p / $maxPoor);					# Will make an intensity gradient from black to gray
			$linkCols[$index]	= $gph->colorAllocate($in,$in,0);	# Gradient from black to yellow
			$unlinkCols[$index]	= $gph->colorAllocate($in,0,$in);	# Gradient from black to purple
			$drawCols{$index} = $linkCols[$index];
			$drawCols{-$index} = $unlinkCols[$index];
		}
		$minGray = 3; $maxGray = 50;
		my @greenScale = my @redScale = ();
		for $p ($minGray..$maxGray) {
			my $in = int (255 - 128 * ($p - $minGray) / ($maxGray - $minGray));	# Gradient from gray to white
			$greenScale[$p]	= $gph->colorAllocate(0,$in,0);						# Gradient from half green to green
			$redScale[$p]	= $gph->colorAllocate($in,0,0);						# Gradient from half red to red
		}
		for ($p = $maxPoor + $stepFactor; $p <= 100; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$gind = ($p > $maxGray) ? $maxGray : int(0.5 + $p);
			$drawCols{$index}	= $greenScale[$gind];
			$drawCols{-$index}	= $redScale[$gind];
		}
	} else {
		for ($p = 0; $p <= 1.01; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$drawCols{$index}	= $black;
			$drawCols{-$index}	= $black;
		}
		$green1	= $gph->colorAllocate(0,51,0);
		$blue1	= $gph->colorAllocate(0,0,51);
		for ($p = 1; $p <= 2.01; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$drawCols{$index}	= $green1;
			$drawCols{-$index}	= $blue1;
		}
		$green2	= $gph->colorAllocate(0,153,0);
		$blue2	= $gph->colorAllocate(0,0,153);
		for ($p = 2; $p <= 3.01; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$drawCols{$index}	= $green2;
			$drawCols{-$index}	= $blue2;
		}
		$green3	= $gph->colorAllocate(0,204,0);
		$blue3	= $gph->colorAllocate(0,0,204);
		for ($p = 3; $p <= 4.01; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$drawCols{$index}	= $green3;
			$drawCols{-$index}	= $blue3;
		}
		$green4	= $gph->colorAllocate(0,255,0);
		$blue4	= $gph->colorAllocate(0,0,255);
		for ($p = 4; $p <= 150; $p += $stepFactor) {
			my $index = int(0.5  + $p / $stepFactor);
			$drawCols{$index}	= $green4;
			$drawCols{-$index}	= $blue4;
		}
	}
	
	my %markPos		= ();
	for my $i (0..$#{$toTest}) {
		my $n			= ${$toTest}[$i];
		$markPos{ $n }	= $i;
		my $x			= $i * $hT + $Border;
		my $y			= $height - $Border + 2 + $wT * length($n);
		$gph->stringUp(gdTinyFont,$x,$Border-2	,$n,$black);
		$gph->stringUp(gdTinyFont,$x,$y			,$n,$black);
		$y 				= $x;
		$x 				= $Border - 2 - $wT * length($n);
		$gph->string(gdTinyFont,$x				,$y,$n,$black);
		$gph->string(gdTinyFont,$width-$Border+2,$y,$n,$black);
	}


	# # # # # # # # # #	Cycle down the rows			= $a -> $aMarker
	for my $a (0..($#{$toTest} - 1)) {				# Row Number
		my $aMarker = ${$toTest}[$a];
		my $yl = $a * $hT + $Border;
		my $yr = $yl + $hT - 1;
		# # # # # # # #	Cycle across the columns	= $b -> $bMarker
		for my $b (($a+1)..$#{$toTest}) {			# Column number
			my $bMarker = ${$toTest}[$b];
			next if ($aMarker eq $bMarker);
			my $xl = $b * $hT + $Border;
			my $xr = $xl + $hT - 1;
			my @Output = ();
			my $totTwo = $myVectors{$aMarker} . $myVectors{$bMarker};
			my @tmpA = ($totTwo =~ /(2)/g);						# Count the combined number of twos
			$totTwo = $#tmpA + 1;
			if ($totTwo > $maxAmbig) {						# Too many Ambiguous
				$gph->filledRectangle($xl+2,$yl+2,$xr-2,$yr-2,$pink);
				$gph->filledRectangle($yl+2,$xl+2,$yr-2,$xr-2,$pink);
				next;
			}
			$tempTime = time;
			my ($theLOD, $devLOD) = split "\t", Probability::twoVector($myVectors{$aMarker}, $myVectors{$bMarker});
			$fpMainLoopTime += time - $tempTime;
			my $index = int($theLOD / $stepFactor);
			my $devIn = int($devLOD / $stepFactor);
			$gph->filledRectangle($xl,$yl,$xr,$yr,$drawCols{$index});
			$gph->filledRectangle($yl,$xl,$yr,$xr,$drawCols{$index});
			$gph->filledRectangle($xl+2,$yl+2,$xr-2,$yr-2,$drawCols{$devIn});
			$gph->filledRectangle($yl+2,$xl+2,$yr-2,$xr-2,$drawCols{$devIn});
		}
	}
	open (FILE, ">$DataLoc/$saveFile.$$.gif") or die "Failure to write to $saveFile.gif: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$DataLoc/$saveFile.$$.gif"><BR>\n);
	&PREDICTTERM($numMrk, $saveFile);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub GRIDPREDICT {
	my $toTest		= shift;	# [0] - ref to array of known vectors that will be tested
	my $saveFile	= shift;	# [1] - filename where data should be saved
	my $totMark		= $#{$toTest} +1;
	my $numMrk = ($totMark - 1) * $totMark / 2;
	&PREDICTINIT($numMrk, $saveFile);
	&GRIDPLOT($toTest);
	my $gridLimit		= 100;
	my $tagTime			= 60 * 5;
	my %copyNumber		= %copyEvnts = ();
	if ($totMark > $gridLimit) {
		print "<FONT COLOR=darkgreen>A table is only drawn for less than $gridLimit markers. \n";
		print "Check the save file to retrieve your data.</FONT><BR>\n";
		print "<PRE>";
	} elsif ($linkStats) {
		print "<FONT COLOR=darkgreen>You have requested statistics only, a marker table will not be drawn.</FONT><BR>\n";
		print "<PRE>";
	} else {
		print "<PRE>";
		print "<TABLE BORDER=1><TR><TD>";
		for my $b (1..$#{$toTest}) {			# Column number
			print "<TH>${$toTest}[$b]";
		}
		print "<TD><TH>Copies\n";
	}
	# # # # # # # # # #	Cycle down the rows			= $a -> $aMarker
	for my $a (0..($#{$toTest} - 1)) {
		my $aMarker = ${$toTest}[$a];
		unless ($totMark > $gridLimit || $linkStats) {
			print "<TR><TH>$aMarker\n";	
			print "<TD COLSPAN=$a>" if ($a > 0);	# Pad the appropriate number of columns as you move down the rows
		} elsif (time - $lastTime > $tagTime) {		# For longer runs add time tags to help track the data
			open (FILE, ">>$DataLoc/$saveFile") or die "PREDICTEND: Failure to append time mark in $saveFile: $!\n";
			my $tp = time - $pTime;
			&CALCTIME(\$tp, \$units);
			print FILE "Time Point\t$pTime\t= $tp $units\n";
			close FILE;	
			print "$tp $units\n";
			$lastTime = time;
		}
		# # # # # # # #	Cycle across the columns	= $b -> $bMarker
		for my $b (($a+1)..$#{$toTest}) {			# Column number
			unless ($totMark > $gridLimit || $linkStats) {
				print "<TD>";
			}
			my $bMarker = ${$toTest}[$b];
			my @Output = ();
			my $totTwo = $myVectors{$aMarker} . $myVectors{$bMarker};
			my @tmpA = ($totTwo =~ /(2)/g);						# Count the combined number of twos
			$totTwo = $#tmpA + 1;
			if ($totTwo > $maxAmbig) {							# Too many Ambiguous
				print "<FONT COLOR=purple>Combined ambiguous results of $totTwo.</FONT>\n" unless ($totMark > $gridLimit || $linkStats);
				next;
			}
			&FASTPREDICT($aMarker,$bMarker, \@Output, $saveFile);
			($bestScore) = split "\t", $Output[0];
			my $avDist = 0; my $avNum = 0;
			foreach $j (@Output) {
				my @d = split "!!", $j;
				if ($d[0] - $bestScore <= $lodCutoff) {
					my @tmpA = ($d[4] =~ /(A)/g);				# Count the number of As
					$copyNumber{$aMarker}[$#tmpA+1]++;
					$copyEvnts{$aMarker}++;
					my @tmpB = ($d[4] =~ /(B)/g);				# Count the number of Bs
					$copyNumber{$bMarker}[$#tmpB+1]++;
					$copyEvnts{$bMarker}++;
					my $dist = $d[2];
					$dist = 1000 if ($dist =~ /unl/i);
					$avDist += $dist; $avNum++;
				}
			}
			&GRIDCELL(\@Output) unless ($totMark > $gridLimit || $linkStats);
#			$avDist = $avDist/ $avNum;
#			if ($avDist > $maxLinkage) {
#				$avDist = $#col;
#			} else {
#				$avDist = $avDist/ $colDiv;
#			}
#			my $x = $Border + $markPos{ $bMarker } * $hT;
#			my $y = $Border + $markPos{ $aMarker } * $hT;
#			$gph->filledRectangle($x,$y,$x+$hT-2,$y+$hT-2,$col[ $avDist ]);
#			$gph->filledRectangle($y,$x,$y+$hT-2,$x+$hT-2,$col[ $avDist ]);
		}
		unless ($totMark > $gridLimit || $linkStats) {
			print "<TH>$aMarker\n";	
			if ($a == 0) {
				print "<TD><PRE>";					# The first marker does not have a column, so copy number must be tacked on to the end
				for my $cn (0..$#{$copyNumber{$aMarker}}) {
					next if ($copyNumber{$aMarker}[$cn] <= 0);
					print $cn + 1;
					my $p = int (100 * $copyNumber{$aMarker}[$cn] / $copyEvnts{$aMarker});
					print ": $p%\n";
				}
			}
		}
	}
	unless ($totMark > $gridLimit || $linkStats) {
		print "<TR><TD>";
		for my $b (1..$#{$toTest}) {			# Column number
			print "<TH>${$toTest}[$b]";
		}
		print "<TR><TH>Copies:\n";
	}
	for my $b (1..$#{$toTest}) {	# Column number
		my $mrk = ${$toTest}[$b];
		print "<TD><PRE>" unless ($totMark > $gridLimit || $linkStats);
		for my $cn (0..$#{$copyNumber{$mrk}}) {
			$copyNumber{$mrk}[$cn] = int (0.5 + 100 * $copyNumber{$mrk}[$cn] / $copyEvnts{$mrk});
			next if ($copyNumber{$mrk}[$cn] <= 0);
			print $cn unless ($totMark > $gridLimit || $linkStats);
			print ": $copyNumber{$mrk}[$cn]%\n" unless ($totMark > $gridLimit || $linkStats);
		}
	}
	open (FILE, ">>$DataLoc/$saveFile") or die "GRIDPREDICT: Failure to append to $saveFile: $!\n";
	for my $b (0..$#{$toTest}) {	# Column number
		my $mrk = ${$toTest}[$b];
		print FILE "Copy\t$mrk";
		for my $cn (0..$#{$copyNumber{$mrk}}) {
			print FILE "\t$copyNumber{$mrk}[$cn]";
		}
		print FILE "\n";
	}
	close FILE;	
	print "</TABLE>" unless ($totMark > $gridLimit || $linkStats);
	print "</PRE>";
	&PREDICTTERM($numMrk, $saveFile);
	if ($linkStats) {
		&STATOUT($saveFile);
	} else {
		&PREDICTEND($saveFile, $timeOut);
	}
	
	open (FILE, ">$DataLoc/$saveFile.$$.gif") or die "Failure to write to $saveFile.gif: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$DataLoc/$saveFile.$$.gif"><BR>\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Print out the results of a comparison in the cell of a grid search
sub GRIDCELL {
	my $outRef		= shift;		# [0] - reference to array holding data
	print "<PRE><FONT COLOR=green>";
	my $w = 0;
	foreach my $n (@{$outRef}) {							# Find maximum length of model text (used for formating output
		my @d = split "!!", $n;
		$w = length($d[4]) if (length($d[4]) > $w);
	}
	my $unlCol = "black";
	for my $j (0..$#{$outRef}) {
		my @d = split "!!", ${$outRef}[$j];
		if ($d[0] - $bestScore > $lodCutoff) {
			print "<FONT COLOR=red>";
			$unlCol = "gray";
		} else {
			print "<FONT COLOR=blue>" if ($j == 1);
		}
		if ($d[2] =~ /unl/i) {
			printf("<FONT COLOR=$unlCol>%${w}.${w}s %5.1f Unlnk</FONT>", $d[4], $d[0]);
		} else {
			printf("%${w}.${w}s %5.1f %3.fkb", $d[4], $d[0], $d[2]);
		}
		print "\n" if ($j != $#{$outRef});
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub NEIGHBORPREDICT {
	my $toTest		= shift;	# [0] - ref to array of known vectors that will be tested
	my $saveFile	= shift;	# [1] - filename where data should be saved
	&PREDICTINIT($#{$toTest}, $saveFile);
	if ($linkStats) {
		print "<PRE>";
	} else {
		print "<TABLE BORDER=1><TR><TH>Marker Pair<TD><PRE>";
		printf("%30.30s  %12.12s at %10.10s kb", "Model", "log Probability", "Distance");
		print "<TH>Status";
	}
#	print "</TABLE>";
	undef(@copyProb);
	for $c (1..$maxCopyNum) {
		$probPos = 1 - Probability::simplePower((1-$RF/100), $c);
		for $k (0..$numHybrids) {
			my $p = Probability::simplePower($probPos, $k) * Probability::simplePower(1 - $probPos, $numHybrids - $k);
			$p *= &NCHOOSEK($numHybrids, $k);
#			print "$p ";
			$copyProb[$c][$k] = -log($p) / log(10);
#			print "$copyProb[$c][$k]\t";
		}
		print "<TH>$c";
#		print "<BR>";
	}


	for my $i (1..$#{$toTest}) {
		my $aMarker = ${$toTest}[$i-1];
		my $bMarker = ${$toTest}[$i];
		print "<TR><TD>$aMarker <FONT COLOR=lightblue>vs.</FONT> $bMarker\n<TD>" unless ($linkStats);
		my $totTwo = $myVectors{$aMarker} . $myVectors{$bMarker};
		my @tmpA = ($totTwo =~ /(2)/g);						# Count the combined number of twos
		$totTwo = $#tmpA + 1;
		if ($totTwo > $maxAmbig) {							# Too many Ambiguous
			print "<FONT COLOR=purple>Combined ambiguous results of $totTwo.</FONT>\n" unless ($linkStats);
			next;
		}
		my @Output = ();
		&FASTPREDICT($aMarker,$bMarker, \@Output, $saveFile);
		next if ($linkStats);
		print "<PRE><FONT COLOR=green>";
		my ($bestScore,$d1,$d2,$d3,$bestMod) = split "!!", $Output[0];
#		print ">> $bestMod";
		my $unlCol = "black";
		my $delta = 0;
		my $state = "Unlinked";
		if ($bestMod =~ /\*/) {
			$state = "Linked";
		}
		for my $j (0..$#Output) {
			@d = split "!!", $Output[$j];
			if ($d[0] - $bestScore > $lodCutoff) {
				print "<FONT COLOR=red>";
				$unlCol = "gray";
			} elsif ($j == 1) {
				print "<FONT COLOR=blue>";
			}
			if ($d[4] !~ /\*/) {
				printf("<FONT COLOR=$unlCol>%30.30s  %5.2f � %4.2f %16.16s</FONT>", $d[4], $d[0], $d[1], "Unlinked Model");
				$delta = $d[0] - $bestScore if ($state eq "Linked" && $delta == 0);
			} else {
				printf("%30.30s  %5.2f � %4.2f at %4.0f � %3.0f kb", $d[4], $d[0], $d[1], $d[2], $d[3] );
				$delta = $d[0] - $bestScore if ($state eq "Unlinked" && $delta == 0);
			}
			print "\n" if ($j != $#Output);
		}
		print "<TD>";
		if ($delta && $delta < $lodCutoff) {
			my $col = ($state eq "Linked" ? "LIME" : "PINK");
			printf("<FONT COLOR=$col>%s %.2f\n", $state, $delta);
		} else {
			my $col = ($state eq "Linked" ? "GREEN" : "RED");
			print "<FONT COLOR=$col><B>$state</B>";
			
		}
		
		my $vec = $myVectors{$aMarker};
		my @tmp2 = ($vec =~ /(2)/g);						# Count the number of twos
		my @tmp1 = ($vec =~ /(1)/g);						# Count the number of ones
		my $hits = int( ($#tmp1 + 1) + ($#tmp2 + 1)/2);
		
		my @cop;
		my $min = 9999999;
		for $c (1..$maxCopyNum) {
			$cop[$c] = $copyProb[$c][$hits];
			$min = $cop[$c] if ($min > $cop[$c]);
		}
		for $c (1..$maxCopyNum) {
			if (abs($min - $cop[$c]) < 3) {
				print "<TD BGCOLOR=yellow><FONT COLOR=green><B>";
			} else {
				print "<TD><FONT COLOR=gray>";
			}
			printf("%.2f", $cop[$c]);
		}
		
		print "\n";
	}
	print "</TABLE>\n";
	&PREDICTTERM($#{$toTest}, $saveFile);
	&STATOUT ($saveFile) if ($linkStats);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PREDICTTERM {
	my $numComp		= shift;	# [0] - number of comparisons that will be done
	my $saveFile	= shift;	# [1] - filename where data should be saved
	$pTime = time - $pTime;
	$pTime = 1 if ($pTime < 1);
	my ($rate, $mkSpeed);
	my $MLR = int (1000 * $fpMainLoopTime / $pTime)/10;
	if ($pTime == 0) {
		$rate = "an awful lot of";
	} else {
		$Preferences{modelRate}			= $rate		= $uniqueVariants * $numSubModels / $pTime;
		$rate		= int (0.5 + $rate);
		$Preferences{avgMarkerSpeed}	= $mkSpeed	= $pTime / $numComp;
		$mkSpeed	= int (100 * $mkSpeed) / 100;
		&WRITEPREF;
	}
	&CALCTIME(\$pTime, \$units);
	print "<BR><FONT COLOR=BROWN><i>";
	$timeOut = "$totalVariants ambiguity variants distilled to $uniqueVariants unique patterns. ";
	$timeOut .= "The average speed was $mkSpeed seconds per marker.<BR>\n";
	$timeOut .= "$numSubModels statistical sub-models (in $numModels models) were analyzed in $pTime $units ";
	$timeOut .= "for a rate of $rate submodel-variants per second.<BR>\n";
	$timeOut .= "The main probability loop occupied $MLR % of the time.<BR>\n";
	print "$timeOut</FONT></i>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub FASTPREDICT {
	my $aMarker		= shift;	# [0] - ref to array of known vectors that will be tested
	my $bMarker		= shift;	# [1] - filename where data should be saved
	my $outRef		= shift;	# [2] - Reference to array that will hold data
	my $saveFile	= shift;	# [3] - filename where data should be saved
	$tempTime = time;
	@{$outRef} = split "\t", Probability::twoVector($myVectors{$aMarker}, $myVectors{$bMarker});
	$fpMainLoopTime += time - $tempTime;
	if ($linkStats) {
		&PREDICTSTAT( $outRef);
	} else {
		&PREDICTADD($saveFile, $aMarker, $bMarker, $outRef);
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub WRITEPREF {
	open (FILE, ">$DataLoc/Preferences") or die "Failure to write Preferences file: $!\n";
	foreach my $k (keys %Preferences) {
		print FILE "$k\t$Preferences{$k}\n";
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub READPREF {
	print "<!"; my $checkPrf = `ls $DataLoc/Preferences`;	chomp $checkPrf; print ">";
	unless ($checkPrf eq "") {
		open (FILE, "$DataLoc/Preferences") or die "Failure to write Preferences file: $!\n";
		while (<FILE>) {
			my @t = split "\t";
			$Preferences{$t[0]} = $t[1];
		}
		close FILE;
	}
	$Preferences{avgMarkerSpeed}	= 0.2 unless ($Preferences{avgMarkerSpeed});
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub XYPLOT {
	my $label		= shift;					# [0] - ref to array of text to use as chart label
	my $Yref		= shift;					# [1] - ref to array of (refs of hash of values)
	my $Xref		= shift;					# [2] - ref to array of X values
	my $width		= shift;					# [3] - width
	my $height		= shift;					# [4] - height
	my $fileName	= shift;					# [5] - filename base
	$fileName		.= int (1000 * rand(1));						# Needed to prevent browser caching
	my $right		= $wT*4	+ 5;
	my $left		= $wT*4	+ 5;
	my $bottom		= $hT*4	+ 5;
	my $top			= 5;
	my $tickSize	= 3;
	my $gph			= new GD::Image($width + $left+$right,$height+$bottom+$top);
	my $white		= $gph->colorAllocate(255,255,255);
	my $black		= $gph->colorAllocate(0,0,0);
	my $red			= $gph->colorAllocate(255,0,0);
	my $green		= $gph->colorAllocate(0,128,0);
	my $blue		= $gph->colorAllocate(0,0,255);
	my $orange		= $gph->colorAllocate(255,128,0);
	my $purple		= $gph->colorAllocate(255,0,128);
	my $gray		= $gph->colorAllocate(150,150,150);
	$gph->transparent($white);
	$gph->interlaced('true');
	my @Xdata 		= @{$Xref};
	my $xmin		= $Xdata[0];
	my $xmax		= $Xdata[$#Xdata];
	my @Ydata		= ();
	my @allKeys		= (); my %inKeys = ();
	for my $s (0..$#{$Yref}) {
		foreach my $n (keys %{${$Yref}[$s]}) {
			next if ($#{$Yref} > 0 && $n eq "$NT$NT");				# Don't include double negatives in multiple plots
			unless (exists $inKeys{$n}) {
				push @allKeys, $n;
				$inKeys{$n} = 1;
			}
			for my $e (0..$#{${$Yref}[$s]{$n}}) {
				$Ydata[$s]{$n}[$e] = ${$Yref}[$s]{$n}[$e];
			}
		}
	}
	my %keyCols		= ();
	my $counter		= 0;
	my $dum			= "";
	my @nt			= (1,2,5,10);									# "pleasing" multiples for tick labeling
	my $ymin		= 100;
	my $ymax		= 0;
	my $NT			= $nullToken;
	my @Ykeys		= ("A$NT", "B$NT", "AB", "$NT$NT", "A$NT & B$NT");
	my @cols		= ($blue, $red, $green, $orange, $purple);

	######## Set up Left Y axis ########
	foreach my $n (@Ykeys) {
		$keyCols{$n} = $cols[$counter++];
		next if ($n eq "$NT$NT");
		$dum .= " $n";
		for my $s (0..$#Ydata) {
			next unless (exists $Ydata[$s]{$n});
			for my $i (0..$#{$Ydata[$s]{$n}}) {
				$Ydata[$s]{$n}[$i] *= 100;
				$ymin = $Ydata[$s]{$n}[$i] if ($ymin > $Ydata[$s]{$n}[$i]);
				$ymax = $Ydata[$s]{$n}[$i] if ($ymax < $Ydata[$s]{$n}[$i]);
			}
		}
	}
	my $maxTick		= $height / ($hT+2);							# The maximum number of y-axis values that can legibily be displayed
	my $tickStep	= int (($ymax-$ymin) / $maxTick)+1;				# The spacing that such ticks would have
	my $factor 		= 0.01;			my $z = 0;

	while ($nt[$#nt] * $factor < $tickStep) { $factor *= 10; }		# Find the order of magnitude that encompases the tickstep
	while ($nt[$z]   * $factor < $tickStep) { $z++;}				# Find the appropriate "nicetick" that suits tickstep
	my $stepFactor = $nt[$z] * $factor;
	$ymin = $stepFactor * int ($ymin / $stepFactor);
	my $yRat = $height / ($ymax-$ymin);
	my ($x, $y);
	for (my $k = $ymin; $k <= $ymax; $k += $stepFactor) {
		$y = $height + $top - ($k-$ymin) * $yRat;
		$gph->line($left-$tickSize,$y,$left-1,$y,$black);
		$gph->string(gdTinyFont,$left-$tickSize-1-$wT*length($k),$y-$hT/2,$k,$black);
	}
	$gph->line($left-1,$top,$left-1,$height+$top,$black);
	$axisLabel = "Probability (%)";
	$axisY = $top + ($height+$wT*length($axisLabel.$dum))/2;
	$gph->stringUp(gdTinyFont,0,$axisY,$axisLabel,$black);
	$axisY -= $wT*length($axisLabel);
	foreach my $n (@allKeys) {
		next if ($n eq "$NT$NT");
		$gph->stringUp(gdTinyFont,0,$axisY," $n",$keyCols{$n});
		$thisYmin{$n} = $ymin;
		$thisYrat{$n} = $yRat;
		$axisY -= $wT*length(" $n");
	}

	if ($#Ydata == 0) {
		######## Set up Right Y axis ########
		$ymin		= 100;
		$ymax		= 0;
		for my $s (0..$#Ydata) {
			next unless (exists $Ydata[$s]{"$NT$NT"});
			for my $i (0..$#{$Ydata[$s]{"$NT$NT"}}) {
				$Ydata[$s]{"$NT$NT"}[$i] *= 100;
				$ymin = $Ydata[$s]{"$NT$NT"}[$i] if ($ymin > $Ydata[$s]{"$NT$NT"}[$i]);
				$ymax = $Ydata[$s]{"$NT$NT"}[$i] if ($ymax < $Ydata[$s]{"$NT$NT"}[$i]);
			}
		}
		if ($ymin == $ymax) {
			$ymin--; $ymax++;
		}
		$thisYmin{"$NT$NT"} = $ymin;
		$maxTick		= $height / ($hT+2);							# The maximum number of y-axis values that can legibily be displayed
		$tickStep		= int (($ymax-$ymin) / $maxTick)+1;				# The spacing that such ticks would have
		$factor 		= 0.01;			my $z = 0;

		while ($nt[$#nt] * $factor < $tickStep) { $factor *= 10; }		# Find the order of magnitude that encompases the tickstep
		while ($nt[$z]   * $factor < $tickStep) { $z++;}				# Find the appropriate "nicetick" that suits tickstep
		$stepFactor = $nt[$z] * $factor;
		$ymin = $stepFactor * int ($ymin / $stepFactor);
		$thisYrat{"$NT$NT"} = $height / ($ymax-$ymin);
		my ($x, $y);
		for (my $k = $ymin; $k <= $ymax; $k += $stepFactor) {
			$y = $height + $top - ($k-$ymin) * $thisYrat{"$NT$NT"};
			$gph->line($left+$width+1,$y,$left+$width+1+$tickSize,$y,$black);
			$gph->string(gdTinyFont,$left+$width+1+$tickSize+1,$y-$hT/2,$k,$black);
		}
		$gph->line($left+$width+1,$top,$left+$width+1,$height+$top,$black);
		$axisLabel = "Probability (%)";
		$axisY = $top + ($height+$wT*length($axisLabel." $NT$NT"))/2;
		$gph->stringUp(gdTinyFont,$width+$left+$right-$hT,$axisY,$axisLabel,$black);
		$axisY -= $wT*length($axisLabel);
		$gph->stringUp(gdTinyFont,$width+$left+$right-$hT,$axisY," $NT$NT",$keyCols{"$NT$NT"});
	}
	
	######## Set up X axis ########
	my $xRat;
	my @tags = ('a','b','c','d','e','f','g','h','i','j','k','l','m','n');
	if ($#Xdata == 0) {
		$axisLabel = "Probabilty Invariant";
		$yLabel = $top+$height+$hT;
		$Xdata[1] = 1; $Xdata[0] = 0;
		for my $s (0..$#Ydata) {
			foreach my $n (keys %{$Ydata[$s]}) {
				$Ydata[$s]{$n}[1] = $Ydata[$s]{$n}[0];
			}
		}
		$xRat = $width;
	} else {
		$maxTick		= $width / ($hT+2);								# The maximum number of x-axis values that can legibily be displayed
		$tickStep		= int (($xmax-$xmin) / $maxTick)+1;				# The spacing that such ticks would have
		$factor 		= 1;			$z = 0;
		while ($nt[$#nt] * $factor < $tickStep) { $factor *= 10; }		# Find the order of magnitude that encompases the tickstep
		while ($nt[$z]   * $factor < $tickStep) { $z++;}				# Find the appropriate "nicetick" that suits tickstep
		$stepFactor = $nt[$z] * $factor;
		$xmin = $stepFactor * int ($xmin / $stepFactor);
		$xRat = $width / ($xmax-$xmin);
		for (my $k = $xmin; $k <= $xmax; $k += $nt[$z] * $factor) {
			$x = $left + ($k-$xmin) * $xRat;
			$gph->line($x,$height+$top+1,$x,$height+$top+$tickSize,$black);
			$gph->stringUp(gdTinyFont,$x-($hT/2),$height+$top+$tickSize+1+ $wT*length($k),$k,$black);
		}
		$axisLabel = "Variable Distance (kb)";
		$yLabel = $top+$bottom+$height-$hT;
	}
	$gph->string(gdTinyFont,$left + ($width-$wT*length($axisLabel))/2,$yLabel,$axisLabel,$black);
	$gph->line($left,$height+1+$top,$width+$left,$height+1+$top,$black);

	$counter = 0; my $keyX = $left*2;
	my @plotZone = ();
	foreach my $n (@allKeys) {
#		$gph->string(gdTinyFont,$keyX,$counter,$n,$keyCols{$n});
		for my $s (0..$#Ydata) {
			next unless (exists $Ydata[$s]{$n});
			my $lasty = $height + $top	- ($Ydata[$s]{$n}[0]	- $thisYmin{$n}) * $thisYrat{$n} + $s;
			my $lastx = $left			+ ($Xdata[0]		- $xmin) * $xRat;
			##### MAJOR FILE LOSS AFTER THIS POINT
			for my $i (0..$#{$Ydata[$s]{$n}}) {
				my $y = $height + $top	- ($Ydata[$s]{$n}[$i]	- $thisYmin{$n}) * $thisYrat{$n} + $s;
				my $x = $left			+ ($Xdata[$i]		- $xmin) * $xRat;
				$gph->line($lastx,$lasty,$x,$y,$keyCols{$n});
				if ($i == $#{$Ydata[$s]{$n}} && $#Ydata > 0) {
					my $looping = 1;
					my $tagx = $x + 3;
					my $tagy = $y - $hT;
					while ($looping) {
						$looping = 0;
						for my $z (0..$#plotZone) {
							my $dx = abs($plotZone[$z][0] - $tagx);
							my $dy = abs($plotZone[$z][1] - $tagy);
							if ($dy < $hT && $dx < $wT) {
								$looping = 1;
								$tagx += $wT;
								$tagy -= 2;
								$tagy = 0 if ($tagy < 0);
								last;
							}
						}
					}
					$gph->string(gdTinyFont,$tagx,$tagy,$tags[$s],$keyCols{$n});
					$gph->dashedLine($x,$y,$tagx,$tagy + $hT/2,$keyCols{$n});
					push @plotZone, [($tagx, $tagy)];
				}
				$lastx = $x;
				$lasty = $y;
			}
			if ($#{$Ydata[$s]{$n}} < $#Xdata) {
				$gph->dashedLine($lastx,$lasty,$width+$left,$lasty,$keyCols{$n});
			}

			$counter += $hT;
		}
	}
	if ($#Xdata > 1) {
		my $x = $left			+ ($maxLinkage		- $xmin) * $xRat;
		$gph->dashedLine($x,$top,$x,$top+$height,$black);
	}
	$wL = (gdGiantFont->width); $hL = (gdGiantFont->height); my $lY = 0;
	for my $lb (0..$#{$label}) {
		my $l = ${$label}[$lb];
		$l = "$tags[$lb]:" . $l if ($#Ydata > 0);
		$gph->string(gdGiantFont,$left + ($width-$wL*length($l))/2,$lY,$l,$black);
		$lY += $hL;
	}
	open (FILE, ">$fileName") or die "Failure to write to $fileName: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$fileName"><BR>\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub HISTOGRAM {
	my $gap			= 1;											# Number of pixels between histogram bars
	my $label		= shift;										# [0] - Name of Histogram (Y Axis)
	my $data		= shift;										# [1] - ref to array of values
	my $width		= shift;										# [2] - width
	my $height		= shift;										# [3] - height
	my $fileName	= shift;										# [4] - filename base
	my $leftVal		= 0;
	if ($#{$data} < 0) {
		print qq(<CENTER><FONT COLOR=red>Data array for the histogram <FONT COLOR=blue>"$label"</FONT> is empty.</FONT></CENTER>\n);
		return;
	}
	while (${$data}[$leftVal] == 0 && $leftVal < $#{$data}) { $leftVal++; }
	my $rightVal	= $#{$data};
	while (${$data}[$rightVal] == 0 && $rightVal > $leftVal) { $rightVal--; }
	my $xRat		= int ($width / ($rightVal-$leftVal+1)) + 1;
	$width			= $xRat * ($rightVal-$leftVal+1);				# We'll make the width as close to the request as possible, but want it to be integral number of pixels
	my $barWidth	= $xRat - $gap;
	$fileName .= int (1000 * rand(1));								# Needed to prevent browser caching
	my $right		= 5;
	my $left		= $wT*4	+ 5;
	my $bottom		= $hT*3	+ 5;
	my $top			= 5;
	
	my $tickSize	= 3;
	my $gph			= new GD::Image($width + $left+$right,$height+$bottom+$top);
	my $white		= $gph->colorAllocate(255,255,255);
	my $black		= $gph->colorAllocate(0,0,0);
	my $red			= $gph->colorAllocate(255,0,0);
	my $green		= $gph->colorAllocate(0,128,0);
	my $blue		= $gph->colorAllocate(0,0,255);
	my $orange		= $gph->colorAllocate(255,128,0);
	$gph->transparent($white);
	$gph->interlaced('true');
	
	my $tot = 0; my $max = 0; my $avg = 0;
	for my $n (0..$#{$data}) {
		$tot += ${$data}[$n];										# The total number of events in the data set
		$max = ${$data}[$n] if (${$data}[$n] > $max);				# Keep track of the maximum value for graphical scaling
		$avg += $n * ${$data}[$n]
	}
	$avg = int (100 * $avg/$tot) / 100;
	
	my $maxPer		= 100 * $max / $tot;							# The percentage of hits in the highest category
	my $yRat		= $height / $maxPer;							# Y scaling ratio - multiply by percenatage value to get pixel offset
	my @nt			= (1,2,5,10);									# "pleasing" multiples for tick labeling
	######## Set up Y axis ########
	my $maxTick		= $height / ($hT+2);							# The maximum number of y-axis values that can legibily be displayed
	my $tickStep	= int ($maxPer / $maxTick)+1;					# The spacing that such ticks would have
	my $factor 		= 0.01;			my $z = 0;

	while ($nt[$#nt] * $factor < $tickStep) { $factor *= 10; }		# Find the order of magnitude that encompases the tickstep
	while ($nt[$z]   * $factor < $tickStep) { $z++;}				# Find the appropriate "nicetick" that suits tickstep
	my ($x, $y);
	for (my $k = 0; $k <= $maxPer; $k += $nt[$z] * $factor) {
		$y = $height + $top - $k * $yRat;
		$gph->line($left-$tickSize,$y,$left-1,$y,$black);
		$gph->string(gdTinyFont,$left-$tickSize-1-$wT*length($k),$y-$hT/2,$k,$black);
	}
	$gph->line($left-1,$top,$left-1,$height+$top,$black);
	######## Set up X axis ########
	$maxTick		= $width / ($hT+2);								# The maximum number of x-axis values that can legibily be displayed
	$tickStep		= int (($rightVal-$leftVal+1) / $maxTick)+1;		# The spacing that such ticks would have
	$factor 		= 1;			$z = 0;
	while ($nt[$#nt] * $factor < $tickStep) { $factor *= 10; }		# Find the order of magnitude that encompases the tickstep
	while ($nt[$z]   * $factor < $tickStep) { $z++;}				# Find the appropriate "nicetick" that suits tickstep
	for (my $k = $leftVal; $k <= $rightVal; $k += $nt[$z] * $factor) {
		$x = $left + ($k-$leftVal) * $xRat + int($barWidth/2);
		$gph->line($x,$height+$top+1,$x,$height+$top+$tickSize,$black);
		$gph->stringUp(gdTinyFont,$x-($hT/2),$height+$top+$tickSize+1+ $wT*length($k),$k,$black);
	}
	$gph->line($left,$height+1+$top,$width+$left,$height+1+$top,$black);
	
	######## Draw the data ########
	for my $k ($leftVal..$rightVal) {
		$x = $left + ($k-$leftVal) * $xRat;
		$y = $height+$top - ($yRat * ${$data}[$k] * 100 / $tot);
		$gph->filledRectangle($x,$y,$x+$barWidth-1,$height+$top,$orange);
	}
	my $labelY = 0; $labelStep = 15;
	my @labelLines = split '\n', $label;
	foreach $n (@labelLines) {
		if ($labelY > 0) {
			$gph->string(gdTinyFont,$left+1,$labelY,$n,$green);
		} else {
			$gph->string(gdSmallFont,$left+1,$labelY,$n,$green);
		}
		$labelY += $labelStep;
	}
	$gph->string(gdTinyFont,$left+1,$labelY,"Average = $avg",$blue);	

	open (FILE, ">$fileName") or die "Failure to write to $fileName: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$fileName"><BR>\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - ref to results hash
# [1] - ref to array of markers to be tested
sub DRAWVECTORS {
	my $results = shift;
	my $marks = shift; my $spaces = "          ";
	print "<TABLE><TR><TD>Markers<TD>\n<PRE>";
	for (my $i = 1; $i < $numHybrids; $i += 10) {
		print "$i";
		print substr ($spaces, 0, 10-length($i));
	}
	print "</TD>\n";
	foreach my $n (@{$marks}) {
		print "<TR><TD>$n<TD>\n<PRE>";
		for my $i (1..$numHybrids) {
			print "${$results}[$i]{$n}";
		}
		print "</TD>\n";
	}
	print "</TABLE>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub STARTHTML {
	print "<html><TITLE>Radiation Hybrid Simulator</TITLE><BODY BGCOLOR=White>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PRINTOPTIONS {
	print "<TABLE BORDER=1>";
	print "<TR><TD ALIGN=center COLSPAN=4 BGCOLOR=GOLD><B><FONT SIZE=+2>Radiation Hybrid Simulator</FONT></B><BR>";
	print "<i><FONT COLOR=olive>Mouse-over any field for help on the functioning of each setting</i>";	
	my $colWidth = 200; my $doubleCol = 2 * $colWidth;
	
	print "<TR><TD VALIGN=TOP COLSPAN=2 WIDTH=$doubleCol>";
	print $query->start_multipart_form(
								-action=>"SimulateRH.cgi",
								-method=>POST,
								-name=>"MainForm");
	print "<FONT COLOR=blue><B>Model Statistics Generation</B></FONT><BR>Panel: ";
	my @list = @theChr;
	push @list, "All";
	print $query->popup_menu(	-name=>'thePanel',
								-title=>"Select the virtual chromosome / marker set that you wish to use for the simulation.",
								-values=>\@list,
								-default=>'$thePanel'), "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;\n";
	print $query->submit(		-value=>'Execute',
								-title=>"Run the program with the settings and files indicated above.",
								),"<BR>\n";
	print "Scan from ";
	print $query->textfield(	-name=>'scanLow',
								-title=>"These settings allow you to automatically calculate statistics for multiple spacings between the TWO CLOSEST MARKERS in your chromosome. This value represents the smallest spacing.",
								-default=>$scanLow,
								-size=>5),"kb to \n";
	print $query->textfield(	-name=>'scanHigh',
								-title=>"This is the highest separation that the two neighbors will reach.",
								-default=>$scanHigh,
								-size=>5),"kb in \n";
	print $query->textfield(	-name=>'scanStep',
								-title=>"The number of distances to simulate for each model.",
								-default=>$scanStep,
								-size=>5),"steps.\n";
	print $query->checkbox(-name=>'logDist',
								-title=>"If checked, the distances will be distributed logarithmically - otherwise they will be linearly spaced.",
								-checked=>$logDist,
								-label=>'',
								-value=>'1'), "Log Spacing<BR>\n";
	print "<TD VALIGN=TOP WIDTH=$colWidth>\n\n";
	print "kb / cR =";
	print $query->textfield(	-name=>'kbPercR',
								-title=>"The number of kilobases per centiRay in the panel.",
								-default=>$kbPercR,
								-size=>6),"<BR>\n";
	print "Retention Frequency ";
	print $query->textfield(	-name=>'RF',
								-title=>"Value between 0 and 100 indicating the average probability of a marker being retained in a hybrid.",
								-default=>$RF,
								-size=>3),"<BR>\n";
	print "Number of Hybrids ";
	print $query->textfield(	-name=>'numHybrids',
								-title=>"The number of hybrids chosen for the panel.",
								-default=>$numHybrids,
								-size=>3),"\n";
	my $nt	= $#vecOrder	+ 1;
	my $n	= $#Statistics	+ 1;
	print "<! "; my $checkOut = `ls $DataLoc/PredictionOutput`;	chomp $checkOut; print ">";
	print "<TD VALIGN=TOP WIDTH=$colWidth>\n\n";
	print qq(<FONT SIZE=-1><A HREF="SimulateRH.cgi?helpMe=Y" target=_blank>Help / FAQ</A></FONT><BR>\n);
	print qq(<FONT COLOR=orange><A HREF="SimulateRH.cgi?ShowStats=Y&dum=$$" target=_blank>View Statistics</A>\tfor $n simulations.</FONT><BR>\n) if ($n > -1);
	print qq(<FONT COLOR=orange><A HREF="SimulateRH.cgi?ShowVectors=Y&dum=$$" target=_blank>View Vectors</A>\tfor $nt markers.</FONT><BR>\n) if ($nt > -1);
	if ($checkOut) {
		my $theSize = -s "$DataLoc/PredictionOutput";
		$theSize = int (0.5 + 100 * $theSize/(1024*1024))/ 100;
		print qq(<FONT COLOR=orange><A HREF="SimulateRH.cgi?Dump=PredictionOutput&dum=$$" target=_blank>View previous Output</A> $theSize Mb.</FONT>\n);
		if ($theSize > 10) {
			print "<BR><FONT COLOR=red SIZE=-1>Too big to view. Drag link to desktop.</FONT>\n";
		}
	}
	print "<TR><TD WIDTH=$colWidth>";
	print "<FONT COLOR=darkgreen><B>Load File</B></FONT>&nbsp;&nbsp;\n";
	print qq( <FONT COLOR=brown SIZE=-1>Click "Execute" to load</FONT><BR>\n);
	print $query->filefield(	-name=>'uploaded_file',
								-default=>$filename,
								-title=>"A text file containing a set of virtual chromosomes",
								-size=>20), "\n";

	print $query->checkbox(-name=>'noStore',
								-title=>"Use for large files - will not attempt to store the file locally, but will instead just read it into memory",
								-checked=>$noStore,
								-label=>'',
								-value=>'1'), "Parse on-the-fly<BR>\n";
	print qq(<FONT COLOR=orange><A HREF="SimulateRH.cgi?seqCompare=Y&dum=$$">Compare to Sequence</A></FONT>  \n);
	print qq(<FONT COLOR=orange SIZE=-2><A HREF="SimulateRH.cgi?seqCompare=Y&noPlace=Y&dum=$$">Plot Only</A></FONT><BR>\n);
	print "Find linked groups LOD => ";
	print $query->textfield(	-name=>'lodLink',
								-title=>"Minimum LOD value to assign a marker to a linkage group.",
								-default=>$lodLink,
								-size=>2),"\n";
		
	my @LodOpt = ("Jones", "Sim");
	print $query->popup_menu(	-name=>'linkMethod',
								-title=>"The method to use when calculating the LOD scores - either the Jones method, or the Simulator's modeling method.",
								-values=>\@LodOpt,
								-default=>'$linkMethod'), "<BR>\n";
	print "<TH>";
	if ($nt) {
		print "Custom Marker Order:<BR>";
		print $query->textarea(		-name=>'myMarkers',
									-default=>$myMarkers,
									-title=>"You can paste a list of marker names (separated by returns) into this field to analyze just a subset of the vector data.",
									-rows=>6,
									-columns=>25);
	}
	print "<TD COLSPAN=2>";
	if ($nt) {
		my @vals = ("No Prediction", "Neighbors", "All Combinations");
		print"Local Vectors: ";
		print $query->popup_menu(	-name=>'doPredict',
									-title=>"If you wish to compare a set of experimental vectors against the simulated library, choose one of these options for best-fit prediction.",
									-values=>\@vals,
									-default=>'$doPredict'), "&nbsp;&nbsp;";
		print"Grid is: ";
		my @gridOpt = ("Ordered", "Alpha");
		print $query->popup_menu(	-name=>'gridOrder',
									-title=>"The order of markers along the grid axes - as found in the file, or alphabetically.",
									-values=>\@gridOpt,
									-default=>'$gridOrder'), "<BR>\n";
		print "<B>Discard</B>: \n";
		print "LOD >";
		print $query->textfield(	-name=>'lodCutoff',
									-title=>"Maximum LOD value (compared to best probability) to use when considering potential next-best model matches.",
									-default=>$lodCutoff,
									-size=>2),"&nbsp;&nbsp;&nbsp;\n";
		print "Ambiguous >";
		print $query->textfield(	-name=>'maxAmbig',
									-title=>"Maximum number of ambiguous assays that can be found in a pairwise comparison.",
									-default=>$maxAmbig,
									-size=>2),"&nbsp;&nbsp;&nbsp;\n";
		print "Distance >";
		print $query->textfield(	-name=>'maxLinkage',
									-title=>"Maximum distance that should be considered for allowing linkage to a model",
									-default=>$maxLinkage,
									-size=>4),"kb.<BR>\n";
		print "Show ";
		print $query->textfield(	-name=>'keepNextBest',
									-title=>"The number of next-best matches to show that have exceeded the LOD limit specified",
									-default=>$keepNextBest,
									-size=>2)," next-best matches.&nbsp;&nbsp;&nbsp;\n";
		print $query->checkbox(-name=>'tossExtreme',
									-title=>"Hits to models at the maximum linkage distance are assumed to be in fact unlinked and are ignored",
									-checked=>$tossExtreme,
									-label=>'',
									-value=>'1'), "Discard extreme matches<BR>\n";
		print $query->checkbox(-name=>'findLinkage',
									-title=>"Calculate the LOD difference between best linked and unlinked models",
									-checked=>$findLinkage,
									-label=>'',
									-value=>'1'), "Calculate Linkage map<BR>\n";
		my $runTime = ($nt-1) * $Preferences{avgMarkerSpeed};
		&CALCTIME(\$runTime, \$units);
		print "<FONT COLOR=magenta>Estimate <U>$runTime $units</u> to calculate neighbors, ";
		$runTime = ($nt-1) * $nt * $Preferences{avgMarkerSpeed}/2;
		&CALCTIME(\$runTime, \$units);
		print "<u>$runTime $units</u> for all combinations.</FONT>\n";
	}
#	print "<TD VALIGN=TOP WIDTH=$colWidth>";
	
	print $query->endform;
	print "</TABLE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DEFINECONSTANTS {
	$HOME = "/usr/local/etc/httpd/cgi-bin/pagelab/charles";
	$Here = "."."$ENV{'REMOTE_ADDR'}";				# Returns the IP address of the *USER'S* Machine (to allow multiple users)
	$DataLoc =	"/tmp/ChasFasta" . $Here;			# Unique directory for the user in the tmp folder
	$WebLoc =	"/tmp/ChasFasta" . $Here;			# The base URL for the host computer
	srand ( time() ^ ($$ + ($$<<15)) );				# Make a "random" seed for random number generation
	unless (-e $DataLoc) {							# Skip this if the directory already exists
		system ("mkdir $DataLoc");
		print "<FONT COLOR=brick>Temporary folder $DataLoc created for data storage.</FONT><BR>\n";
		system("chmod 777 $DataLoc");
	}
	print "<!"; system "rm $DataLoc/*.gif*";	print ">\n";		# Remove old gif files
	print "<!"; system "rm $DataLoc/*.tmp";	print ">\n";
	$panelsPerSec = 7;											# This helps with time estimates
	$hT = (gdTinyFont->height);
	$wT = (gdTinyFont->width);
	@CGIParams = (	"thePanel=",		"kbPercR=5.3",		"RF=12",			"numHybrids=90",
					"numTrials=1",		"falsePositive=0",	"falseNegative=0",	"nullToken=0",
					"scanLow=",			"scanHigh=",		"scanStep=",		"doPredict=",
					"lodCutoff=3",		"maxAmbig=6",		"keepNextBest=1",	"gridOrder=Ordered",
					"maxLinkage=350",	"tossExtreme=1",	"noStore=",			"myMarkers=",
					"linkToken=*",		"unlinkToken=|",	"logDist=1",		"minCopyMajority=75",
					"seqCompare=",		"findLinkage=",		"maxCopyNum=8",		"lodLink=",
					"linkMethod=Jones",	
				);

	$query = new CGI;
	$filename =	$query->param('uploaded_file');
	foreach $i (@CGIParams) {	
		($var, $val) = split "=", $i;
		$$var = $query->param($var) if ($query->param($var) ne "");
		$$var = $val if ($$var eq "");
		$GV{$var} = $$var;
	}
	$lambda = 0.01 / $kbPercR;
	$GV{lambda} = $lambda;
	$DecimalPlaces = 3;
	$DecMult = 1;
	for my $i (1..$DecimalPlaces) { $DecMult *= 10; }
	$findLinkage = 1 if ($seqCompare);
	my @export = ($maxLinkage, $tossExtreme, $lodCutoff, $maxAmbig, $keepNextBest, $numHybrids, $cDebugger, $findLinkage);
	my $k = "d" . ($#export+1);
	Probability::setValues(pack($k,@export));
	&READPREF;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Netscape does not clean up hex representations of special characters (ie leaves things like ":" as "%3A")
sub TIDYBROWSER {
	undef (@tidy);
	@tidy = split "%", $_[0];	
	$cleaned = $tidy[0];
	for $k (1..@tidy) {
		$n = substr ($tidy[$k],0,2);
		$cleaned .= chr(hex($n)) . substr($tidy[$k],2);
	}
	return $cleaned;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub STNDEV {
	my $dataRef	= shift;		# [0] Reference to array of values
	my $avgRef	= shift;		# [1] Ref to value that will hold the average
	my $devRef	= shift;		# [2] Ref to value that will hold standard deviation
	
	$$avgRef = $$devRef = "N/A";
	my $tot = $#{$dataRef} + 1;
	return if ($tot <= 0);
	foreach my $n (@{$dataRef}) {
		$$avgRef += $n;
#		print "$n + ";
	}
	$$avgRef /= $tot;
#	print " = $$avgRef<BR>";
	foreach my $n (@{$dataRef}) {
		$$devRef += ($n - $$avgRef) * ($n - $$avgRef);
	}
	$$devRef	= sqrt($$devRef/$tot);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Takes a user specified file and reads it into memory
sub DYNAMICLOAD {
	%GeneralValues	= ();
	%CopyNum		= ();
	%majorCopy		= ();
	%result			= ();
	my %unl			= ();
	my %wellBehaved	= ();
	my $lt			= $linkToken;
	my $ut			= $unlinkToken;
	my $ti			= time;
	my $fileLines	= my $contentLines = 0;
	my $fileTime	= $elapsedTime = time;
	print "<PRE>";
	while (<$filename>) {
		my @f = split /[\r\n]/;
		foreach my $n (@f) {
			$fileLines++;
			next if ($n eq "");
			next if ($n =~ /time/i);
			if ($n =~ /start/i) {
				$adding = 1;
				next;
			}
			$adding = 0 if ($n =~ /stop/i);
			if ($adding) {
				$contentLines++;
				my @t		= split "\t", $n;
				my $type	= shift @t;
				if (time - $elapsedTime > 60) {
					$elapsedTime = time;
					my $ft = $elapsedTime - $fileTime;
					&CALCTIME(\$ft, \$units);
					printf("%5.1 %7.7s = %6.d lines\n", $ft, $units, $contentLines);
				}
				if ($t[0] eq "" or $t[1] eq "") {
					print qq(<FONT COLOR=red>Fully or partially undefined key/value pair in file: "$n"</FONT>\n);
					next;
				}

				if ($type =~ /val/i) {						# The line contains a simple value
					$GeneralValues{$t[0]} = $t[1];			# Use the first element as the key, and the second as the value
					next;
				}
				if ($type =~ /cop/i) {						# The line contains copy number information
					my $marker			= shift @t;			# Use the first element as the name of the marker
					$CopyNum{$marker}	= @t;				# The remainder of the line designates the copy number estimates
					my $max = 0; my $index = 0;
					for my $i (0..$#t) {
						if ($t[$i] > $max) {
							$max	= $t[$i];
							$index	= $i;
						}
					}
					if ($max > $minCopyMajority) {
						$majorCopy{$marker} = $index;
					}
					next;
				}
				my $markerA		= $type;
				my $markerB		= shift @t;
				
				my @first	= split '!!', $t[0];
				my $priMod	= $first[$#first];
				my $priProb	= $first[0];
				my $priDist	= $first[2];
				
				for my $i (1..$#t) {
					my @d = split '!!', $t[$i];
					my $thisMod		= $d[$#d];
					push @{$LODs{$priMod}{$thisMod}}, $d[0] - $priProb;
					push @{$Dist{$priMod}{$thisMod}}, $d[2] - $priDist;
				}
				
				
			}
		}
	}
	print qq(</PRE><FONT COLOR=green>File "$filename" has $contentLines lines of content out of $fileLines lines total length.</FONT><BR>\n);
	if ($contentLines == 0) {
		&BROWSERERR;
		return;
	}
	my $avg, $dev;
	
	foreach my $a (keys %LODs) {
		foreach my $b (keys %{$LODs{$a}}) {
			next if ($a eq $b);
			&STNDEV( \@{$LODs{$a}{$b}}, \$avg, \$dev);
			my $col = "";
			if 		($avg >= 3) {
				$col = "<FONT COLOR=green>";
			} elsif ($avg >= 2) {
				$col = "<FONT COLOR=yellow>";
			}
			print "$col$a\t$b";
			if ($avg !~ "N") {
				printf("\t%.2f\t�\t%.2f", $avg, $dev);
			}
			unless ($Invariant{$a} || $Invariant{$b}) {
				&STNDEV( \@{$Dist{$a}{$b}}, \$avg, \$dev);
				printf("\t%.2f\t�\t%.2f", $avg, $dev);
			}
			print "</FONT>" if ($col ne "");
			print "<BR>\n";
		}
	}
	
	exit;
	
	
	my $loadTime = time - $ti;
	$ti = time;
	my @allMarkers = sort {$a cmp $b} keys %CopyNum;
	my $numMarkers = $#allMarkers + 1;
	my %trimmed = ();
	my %working = ();
	foreach my $mA (keys %result) {
		foreach my $mB (keys %{$result{$mA}}) {
			my @t = @{$result{$mA}{$mB}};						# @t will hold the models that are considered
			my @s = ();
			my @first = split '!!', $t[0];
			foreach my $e (@t) {								# Find the surviving models that are within LOD cutoff
				my @d = split '!!', $e;
				unless ($d[0] - $first[0] > $lodCutoff) {		# Compare each model to the first one (highest prob)
					push @s, $e;
				}
			}
			@t = ();
			if ($majorCopy{$mA} > 0 || $majorCopy{$mB} > 0) {	# At least one of the markers has had a single copy number determined
				foreach my $e (@s) {
					my @d = split '!!', $e;
					if ($majorCopy{$mA} > 0) {
						my @tmpA = ($d[4] =~ /(A)/g);			# Count the number of As in the model
						if ($#tmpA+1 != $majorCopy{$mA}) {		# Disregard any models that don't have the appropriate number of markers
							$trimmed{$mA}++;
							next;
						}
					}
					if ($majorCopy{$mB} > 0) {
						my @tmpB = ($d[4] =~ /(B)/g);			# Count the number of Bs in the model
						if ($#tmpB+1 != $majorCopy{$mB}) {		# Disregard any models that don't have the appropriate number of markers
							$trimmed{$mB}++;
							next;
						}
					}
					push @t, $e;								# Models that survive copy number prediction are put into @t
				}
			} else {
				@t = @s;
			}
			if ($#t == 0) {
				push @{$wellBehaved{$first[$#first]}}, "$mA\t$mB";
			}
			@{$working{$mA}{$mB}} = @{$working{$mB}{$mA}} = @t;
		}
	}
	
	
	print "<TABLE><TR><TH COLSPAN=2>General Values\n";
	foreach my $n (keys %GeneralValues) {
		print "<TR><TD>$n<TD>$GeneralValues{$n}\n";
	}
	print "</TABLE>\n";
	
	print "<TABLE><TR><TH COLSPAN=3>Copy Number Model Elimination\n";
	foreach my $k (keys %trimmed) {
		print "<TR><TD>$k<TD>$majorCopy{$k} Copies<TD>$trimmed{$k} models killed\n";
	}
	print "</TABLE>\n";
	
	my @easiest	= ("A${lt}B", "A${lt}B${ut}A${lt}B", "A${lt}B${ut}A${lt}B{ut}A${lt}B");
	my %groups	= my %linked	= ();

	print "<PRE>";
	my @s = ();
	foreach my $e (keys %wellBehaved) {
		my $t = $#{$wellBehaved{$e}} + 1;
		push @s, sprintf ("%6.6s\t%s", $t, $e);
	}
	@s = sort {$b cmp $a} @s;
	foreach my $e (@s) {
		my @t = split "\t", $e;
		printf ("%28.28s = %5.0f\tmembers\n", $t[1], $t[0]);
	}
	print "</PRE>";
	
	foreach my $eas (@easiest) {
		foreach my $e (@{$wellBehaved{$eas}}) {
			my ($mA, $mB);
			($mA, $mB) = split "\t", $e;
			my $lA = $linked{$mA};
			my $lB = $linked{$mB};
			if		($lA ne "" && $lB ne "") {					# both already in a group
				if	($lA != $lB) {								# They are in differnt groups, combine the two groups
					foreach my $tr (@{$groups{$eas}[$lB]}) {
						push @{$groups{$eas}[$lA]}, $tr;
						$linked{$tr} = $lA;
					}
					@{$groups{$eas}[$lB]} = ();
				}
			} elsif	($lA ne "") {
				push @{$groups{$eas}[$lA]}, $mB;
				$linked{$mB} = $lA;
			} elsif	($lB ne "") {
				push @{$groups{$eas}[$lB]}, $mA;
				$linked{$mA} = $lB;
			} else {
				push @{$groups{$eas}}, [($mA, $mB)];
				$linked{$mA} = $linked{$mB} = $#{$groups{$eas}};
			}
		}
	}
	
	my %grouped = ();
	foreach my $n (keys %groups) {
		print "<H3>High quality linkage groups of type $n</H3>\n<PRE>";
		my $count = 0;
		for my $g (0..$#{$groups{$n}}) {
			next if ($#{$groups{$n}[$g]} < 0);
			$count++;
			printf ("<B>%3.3s:</B> ", $count);
			foreach my $m  (@{$groups{$n}[$g]}) {
				printf ("%10.10s ", $m);
				$grouped{$m} = "Y";
			}
			print "\n";
		}
		print "</PRE>";
	}
	my @t = ();
	foreach my $k (@allMarkers) {
		push @t, $k unless ($grouped{$k});
	}
	my $tn = $#t + 1;
	print qq(<FONT COLOR=blue><B>The following $tn markers (out of $numMarkers) were not assigned to any groups:</B>\n);
	print "<PRE>"; my $cc = 1;
	foreach my $n (@t) {
		printf ("%15.15s ", $n);
		print "\n" unless ($cc % 6); $cc++;
	}
	print "</FONT></PRE>";
	my $analysisTime = time - $ti;
	&CALCTIME(\$loadTime, \$units);
	print "<FONT COLOR=brown>$loadTime $units to load and parse file of $fileLines lines.<BR>\n";
	&CALCTIME(\$analysisTime, \$units);
	print "$analysisTime $units to analyze data.<BR></FONT>\n";
	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Takes a user specified file and transfers it to the temp folder
sub EXTERNALLOAD {
	my $toTransfer = ""; my $adding = 0;
	my $limit = 100;
	while (<$filename>) {
		my @t = split /[\r\n]/;
#		print "$#t<BR>";
		foreach my $n (@t) {
#			print "$adding $n<BR>";
			if ($n =~ /start/i) {
				$adding = 1;
				next;
			}
			$adding = 0 if ($n =~ /stop/i);
			$toTransfer .= "$n\n" if ($adding);
			if ($tLength > $limit) {
				print "<FONT COLOR=red><B>File transfer of $tLength kb aborted. Limit of $limit kb on imported file.</FONT></B></BR>\n";
				return;
			}
		}
	}
	my $tLength = int(length($toTransfer) / 1024);
	if ($toTransfer =~ /prob/i) {
		&SAVEFILE(\$toTransfer, "LocalStatistics", "a Hybrid Statistics file");
		$newStats = "Y";
	} elsif ($toTransfer =~ /vector/i) {
		&SAVEFILE(\$toTransfer, "LocalVectors", "a Vector Data file");
		$newVec = "Y";
	} elsif ($toTransfer =~ /model def/i) {
		&SAVEFILE(\$toTransfer, "LocalChromosomes", "a Chromosome Data file");
		$newChr = "Y";
	} elsif ($tLength < 1) {
		&BROWSERERR;
	} else {
		print qq(<FONT COLOR=red>Sorry! I could not interpret what form of data your file contains! Please remember to:</FONT><BR>\n);
		print qq(<UL><LI>Start all files with <FONT COLOR=blue>START</FONT> and end them with <FONT COLOR=blue>STOP</FONT>\n);
		print qq(<LI>Include the tag line "Vector File" in your real-data vector files\n);
		print qq(<LI>Do not use instances of the string fragment "prob" in any files <i>other than</i> statistics files.\n);
		print qq(<LI>Do not use instances of the string fragment "chr" in any files <i>other than</i> chromosome model files.\n);
		print "</UL>\n";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub BROWSERERR {
	print qq(<FONT COLOR=red>No data detected in file! Plese remember to:</FONT><BR>\n);
	print qq(<UL><LI>Start all files with <FONT COLOR=blue>START</FONT> and end them with <FONT COLOR=blue>STOP</FONT>\n);
	print qq(<LI>Macintosh users of internet explorer: IE improperly reads part of the resource fork in addition to the data fork of the file \n);
	print "(this is why I've included START/STOP tags). It also sometimes skips large chunks of data from the front of the file. Try copying all \n";
	print "your data, paste it into a new document, and save it. This often fixes the problem.\n";
	print "</UL>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - Reference to string containing data
# [1] - name of file to store data in
# [2] - type of file transfered (used for user feedback)
sub SAVEFILE {
	open (FILE, ">$DataLoc/$_[1]") or die "Failure to write to $_[1]: $!\n";
	print FILE ${$_[0]};
	close FILE;
	print "<FONT COLOR=darkgreen><B>File $filename transfered to local folder as $_[2].</B></FONT><BR>\n";
	system("chmod 666 $DataLoc/$_[1]");
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADDATA {
	print "<! ";
	my $checkChr = `ls $DataLoc/LocalChromosomes`;	chomp $checkChr;
	my $checkVec = `ls $DataLoc/LocalVectors`;		chomp $checkVec;
	my $checkSts = `ls $DataLoc/LocalStatistics`;	chomp $checkSts;
	print ">";
	unless ($checkChr) {
		print "<FONT COLOR=green>Creating example Chromosome file. You can load your own chromosomes using the dialog box below.</FONT><BR>\n";
		my $demo = "Model Definition\n\n";
		$demo	.= "-\t3000\n*\tVariable\t50\n|\tUnlinked\n\n";
		$demo	.= "A*B\nA*B|B\nA*B|A*B\nA|B\nA|B|B\nA|A|B|B\n";
		open (FILE, ">$DataLoc/LocalChromosomes") or die "Failure to write to demonstration file LocalChromosomes: $!\n";
		print FILE $demo;
		close FILE;
	}
	&LOADCHR;
	&LOADSTS if ($checkSts);
	&LOADVEC if ($checkVec);
	&LOADALIAS("$HOME/Data/GeneAliases");
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADALIAS {
	my $filename =	shift;
	undef(%geneAlias);
	open (FILE, "$filename") or die "Failure to open $filename: $!\n";
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		my ($sts, $gene)		= split "\t";
		next if ($gene eq "BE");
		$geneAlias{$sts} = $gene;
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADVEC {
	undef (%myVectors); undef (@vecOrder);
	open (FILE, "$DataLoc/LocalVectors") or die "Failure to open LocalVectors: $!\n";
	my $defaultSize = 0;
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		my @t		= split "\t";
		next if ($t[0] =~ /vector/i);
		next if ($t[0] eq "");
		$defaultSize = length($t[1]) unless ($defaultSize);
		if ($defaultSize != length($t[1])) {
			printf ("<FONT COLOR=red>Vector for %s is rejected for having %d assays instead of %d.</FONT><BR>\n", $t[0], length($t[1]), $defaultSize);
			next;
		}
		$myVectors{$t[0]} = $t[1];
		push @vecOrder, $t[0];
	}
	close FILE;
	@vecAlpha = sort {$a cmp $b} keys %myVectors;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADSTS {
	undef(@Statistics);
	my %GeneralValues;
	my $modelName		= "";
	my @distances		= ();
	my %probabilities	= ();
	my %tokens			= ();
	$numModels			= 0;
	$numSubModels		= 0;
	open (FILE, "$DataLoc/LocalStatistics") or die "Failure to open LocalStatistics: $!\n";
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		if (/###/) {
			next if ($modelName eq "");						# Initially there will be no assignments made, so ignore the first ### for a series
			my $theRef = chrStat->new( $modelName, \%GeneralValues, \%probabilities, \@distances, \%tokens );
			$numModels++;	$numSubModels += $#distances + 1;
			push @Statistics, $theRef;
			$modelName		= "";
			@distances		= ();
			%probabilities	= ();
			next;										# ... but leave the global values intact - allows the values to be specified once for many entries
		}
		my @t		= split "\t";
		my $type	= shift @t;
		if ($t[0] eq "" or $t[1] eq "") {
			print qq(<FONT COLOR=red>Fully or partially undefined key/value pair in LocalStatistics: "$_"</FONT><BR>\n);
			next;
		}
		SWITCH: {
			if ($type =~ /val/i) {						# The line contains a simple value
				$GeneralValues{$t[0]} = $t[1];			# Use the first element as the key, and the second as the value
				last SWITCH;
			}
			if ($type =~ /mod/i) {						# The line contains the header for a model designation
				$modelName = shift @t;					# Use the first element as the name of the model
				@distances = @t;						# The remainder of the line designates the distances used for calculating the probabilities
				last SWITCH;
			}
			if ($type =~ /prob/i) {						# The line contains a probability array
				my $aName = shift @t;					# Use the first element as the key
				@{$probabilities{$aName}} = @t;			# Set the remainder as the value (array) of the hash
				last SWITCH;
			}
			if ($type =~ /tok/i) {						# The line contains a token assignment
				$tokens{$t[0]} = $t[1];			# Use the first element as the key, and the second as the value
				last SWITCH;
			}
	#		print qq(<FONT COLOR=red>Unrecgonized tag line "$type" in LocalStatistics: $_</FONT><BR>\n);
		}
	}
	close FILE;
	my $order = Probability::modelOrder(void);
	@modOrd = split "\t", $order;
#	foreach my $k (@modOrd) {print qq("$k"<BR>);}
	@alphaModels = ();
	foreach $s (@Statistics) {
		push @alphaModels, $s->{Model};
		$Invariant{$s->{Model}} = "Y" if ($s->{Distance}[0] =~ /inv/i);
	}
	@alphaModels = sort {$a cmp $b} @alphaModels;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - ref to chromosome being analyzed
# [1] - filename where data should be saved
sub SAVEPREP {
	my $saveFile	= shift;
	open (FILE, ">$DataLoc/$saveFile") or die "Failure to write to $saveFile: $!\n";
	print FILE "START\n\n";
	print FILE "Value\tAssays per Panel\t$numTrials\n";
	print FILE "Value\tRetention Frequency\t$RF\n";
	print FILE "Value\tkb per cR\t$kbPercR\n";
	print FILE "Value\tFalse Positives\t$falsePositive\n";
	print FILE "Value\tFalse Negatives\t$falseNegative\n";
	print FILE "\n";
	close FILE;	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - ref to statistics hash
# [1] - ref to chromosome being analyzed
# [2] - name of "MarkerA"
# [3] - name of "MarkerB"
# [4] - filename where data should be saved
sub SAVEDATA {
	my $dataRef		= shift;
	my $distRef		= shift;
	my $chromosome	= shift;
	my $saveFile	= shift;
	my @distances	= @{$distRef};
	my @results		= sort {$a cmp $b} keys %{$dataRef};
	
	open (FILE, ">>$DataLoc/$saveFile") or die "Failure to append to $saveFile: $!\n";
	print FILE "##################################################\n";
	print FILE "Model\t$chromosome->{Name}";
	my $max = $#distances; my $invariant = 1;
	foreach my $n (@results) {
		for my $i (1..$max) {
			if (${$dataRef}{$n}[$i] != ${$dataRef}{$n}[0]) {
				$invariant = 0;
				last;
			}
		}
	}
	if ($invariant) {
		$max = 0;
		@distances = ("Invariant");
	}
	for my $i (0..$max) {
		print FILE "\t$distances[$i]";
	}
	print FILE "\n";
	foreach my $n (@results) {
		print FILE "Probability\t$n";
		for my $i (0..$max) {
			printf FILE ("\t%1.5f",${$dataRef}{$n}[$i]);
		}
		print FILE "\n";
	}
	foreach my $n (keys %{$chromosome->{Token}}) {
		print FILE "Token\t$n\t$chromosome->{Token}{$n}\n";		 
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub CALCTIME {
	my $timeRef = shift;
	my $unitRef	= shift;
	$$unitRef = "second";
	if ($$timeRef > 60) {
		$$unitRef = "minute";
		$$timeRef /= 60;
	}
	if ($$timeRef > 60) {
		$$unitRef = "hour";
		$$timeRef /= 60;
	}
	$$timeRef = int (10 * $$timeRef)/10;
	unless ($$timeRef <= 1) {
		$$unitRef .= "s";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SAVEEND {
	my $saveFile	= shift;
	open (FILE, ">>$DataLoc/$saveFile") or die "SAVEEND: Failure to write to $saveFile: $!\n";
	print FILE "##################################################\n";
	print FILE "\nSTOP";
	close FILE;
	
	print qq(<FONT COLOR=blue>Data saved to file <A HREF="SimulateRH.cgi?Dump=$saveFile&dum=$$" target=_blank>$saveFile</A>.</FONT><BR>\n);
	my $rate;
	if ($theTime == 0) {
		$rate = "an infinite number of";
	} else {
		$rate = int(0.5 + $theOperations/$theTime);
	}
	&CALCTIME(\$theTime, \$units);
	print "<FONT COLOR=brown>$theTime $units to analyze $theOperations models = $rate models/sec.<BR></FONT>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADCHR {
	open (FILE, "$DataLoc/LocalChromosomes") or die "Failure to open LocalChromosomes: $!\n";
	undef (%Chromosomes);
	my %alias	= ();
	my %myToken	= ();
	while (<FILE>) {
		chomp;
		next if ($_ eq "");
		next if (/model def/i);
		my @t = split "\t";
		if ($t[1] ne "") {							# A second element indicates that tokens are being assigned
			if		($t[1] =~ /var/i) {				# Variable position
				$myToken{$t[0]} = "VAR";			# Tokens need only be assigned once, but can be reassigned
				if ($t[2] > 0) {					# The third field is the default value to use for variable markers
					$myToken{VAR}	= $t[2];
				} else {
					$myToken{VAR}	= 50;			# If not specified, just pick 50 for the hell of it
				}
			} elsif	($t[1] =~ /unl/i) {				# Unlinked position
				$myToken{$t[0]} = "UNL";			# The marker is designated as unlinked
				if ($t[2] > 0) {					# The third field is the default value to use for unlinked markers
					$myToken{UNL}	= $t[2];
				} else {
					$myToken{UNL}	= 10 / $lambda;	# Set generic "unlinked" tokens to 10 x average fragment length
				}
			} else {								# Specific distance specified
				$myToken{$t[0]} = &CONVERTTOKB($t[1]);
			}
		} elsif	($t[0] =~ '=') {					# The line is specifying an alias
			$alias{$t[0]} = $t[1];
		} else {									# The line is specifying a chromosome
			my $current = Chromosome->new( $t[0], \%myToken,  \%alias);
			$Chromosomes{ $t[0] } = $current;
			%alias = ();							# Aliases are only for one chromosome
		}
	}
	close FILE;

	@theChr = sort {$a cmp $b} keys %Chromosomes;
	if ($newChr) {						# The chromosome was loaded, give the user a summary of what they entered
		my $counter = 0;
		print "<PRE>";
		foreach my $n (@theChr) {
			$counter++;
			$Chromosomes{$n}->drawChromosome(500,"$DataLoc/chromosome$counter$$.gif", \%GV);
			print "<HR>\n";
		}
		print "</PRE>";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - the value being parsed for genomic units
sub CONVERTTOKB {
	$multiplier = 1;
	$multiplier = 1000 if ($_[0] =~ /Mb/i);
	$multiplier = .001 if ($_[0] =~ /bp/i);
	my @t = split /(m|M|k|K|b|B| )/, $_[0];		# letters and whitespaces
#	print "$_[0] -> $t[0] * $multiplier<BR>";
	return ($t[0] * $multiplier);	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub HELPME {
	my $checkHlp = `ls SRHHelp.html`;	chomp $checkHlp;
	if ($checkHlp) {
		open (FILE, "SRHHelp.html") or die "Failure to open SRHHelp.html: $!\n";
		while (<FILE>) {
			print;
		}
	} else {
		print "<html><TITLE>ERROR - Help file missing</TITLE><BODY BGCOLOR=White>\n";
		print "<H2>Sorry, the file <FONT COLOR=blue>SRHHelp.html</FONT> was not found. ";
		print "It should be in the same directory as the program file.</H2>\n";
	}

exit;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This package manages virtual chromosomes used for radiation hybrid simulation
package Chromosome;
use GD;
sub new {
	my $class	= shift;
	my $name	= shift;
	my $tokRef	= shift;
	my $aliRef	= shift;
	my $self = {};
	$self->{Name} = $name;
	foreach my $n (keys %{$tokRef}) {
		$self->{Token}{$n} = ${$tokRef}{$n};		 
	}
	foreach my $n (keys %{$aliRef}) {
		$self->{Alias}{$n} = ${$aliRef}{$n};		 
	}
	bless $self, $class;
	return $self;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub chrParse {
	my $self	= shift;
	my $markRef	= shift;	# [0] ref to array holding marker names
	my $distRef	= shift;	# [1] ref to array holding distances between markers
	my $totDist	= 0;
	my @toParse	= split '', $self->{Name};
	@{$markRef} = @{$distRef} = ();
	while ($#toParse > -1) {
		push @{$markRef}, (shift @toParse);
		if ($#toParse > -1) {
			my $tok = shift @toParse;
			if ($self->{Token}{$tok} eq "VAR") {			# A variable marker
				push @{$distRef}, "V";
				$totDist += $self->{Token}{VAR};			# Keep track of the total distance
			} elsif ($self->{Token}{$tok} eq "UNL") {		# An unlinked marker
				push @{$distRef}, "U";
				$totDist += $self->{Token}{UNL};			# Keep track of the total distance
			} else {										# Distance pre-defined
				push @{$distRef}, $self->{Token}{$tok};
				$totDist += $self->{Token}{$tok};			# Keep track of the total distance
			}
		}
	}
	push @{$distRef}, 0;									# This keeps the distance array the same size as the marker array, may avoid problems later
	return $totDist;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - imported global values
sub calcProbs {
	my @Distances;
	my @marker = (); my @separation = ();
	my $self		= shift;
	my $gv			= shift;								# [0] Imported global values
	my $distRef		= shift;								# [1] Ref to array of distance values for variable distances
	my $outputRef	= shift;								# [2] Ref to hash that will hold the final variables
	%{$outputRef}	= ();
	my $retFreq		= ${$gv}{RF}/100;
	my $lambda		= ${$gv}{lambda};
	my @binaries	= ();
	my $Break		= " ";
	my $noBrk		= "-";
	my $nullHolder	= ${$gv}{nullToken};
	my $nullString	= $nullHolder; $nullString x= 20;
	my @keyContent	= ($nullHolder . $nullHolder, "A" . $nullHolder, "B" . $nullHolder, "AB");
	@{$binaries[1]} = ($Break,$noBrk);
	my $totDist		= $self->chrParse(\@marker, \@separation);
	if (defined $distRef) {
		@Distances	= @{$distRef};
	} else {
		@Distances	= (0);
	}
	for my $i (2..($#marker+1)) {							# Recursively generate list of all possible breaks and non-breaks
		foreach $n (@{$binaries[$i-1]}) {
			push @{$binaries[$i]}, "$n$Break";
			push @{$binaries[$i]}, "$n$noBrk";
		}
	}
	my @breakList	= sort {$a cmp $b} @{$binaries[$#marker]};
	my %probability	= (); my @variables	= (); my %in;
#	print "<PRE>";
	foreach $n (@breakList) {							# Calculate the probability for each break/no break possibility
		my $prob	= 1; my $breakProb;
		my $output	= "";
		for my $i (0..($#marker-1)) {
			my $sep			= $separation[$i];
			my $breakToken	= substr($n,$i,1);
			if		($sep eq "V") {						# Variable marker, note it on a list for future consideration
				push @variables, $i unless ($in{$i});	# Keep track of variables, but only add them once
				$in{$i} = 1;
			} elsif	($sep eq "U") {						# Forced unlinked, break probability = 100% (eg on separate chromosomes)
				unless ($breakToken eq $Break) {
					$prob = 0;							# Only consider combinations that have a break here
				}
			} else {									# Normal breakage, calculate based on distance separating markers
				$breakProb	= 1 - exp(-$lambda * $sep);
				if ($breakToken eq $Break) {
					$prob *= $breakProb;				# Token specifies break
				} else {
					$prob *= (1-$breakProb);			# No break, prob = 1 - probabilityOfBreak
				}
			}
			$output .= "$marker[$i]" . "$breakToken";
		}
		next if ($prob <= 0);							# Don't bother considering break combinations that have 0% chance of occuring
		$output	.= $marker[$#marker];					# Tack on the name of the last marker for the key
#		print "Breaklist $output  $prob\n";
		$probability{$output} = $prob;					# We'll use as a key a handy visual representation
	}
	@pKeys = sort {$a cmp $b} keys %probability;
	
	######### Calculate the Probabilites #########
	my %Content = ();
	foreach my $k (@pKeys) {								# Look at each break/no break possibility
#		print "ProbCalc $k  $probability{$k}\n";
		my @fragments = split "$Break", $k;
		my $numEvents = $#fragments+1;
		my @eventList= sort {$a cmp $b} @{$binaries[$numEvents]};
		my %subContent = ();
		foreach my $event (@eventList) {					# Use the previously calculate binaries array to look at all possible retained/not retained probabilities
			my $retProb	= 1;
			my @content	= ();
			for my $f (0..$#fragments) {					# Now consider each fragment to see if it is retained
				my $eventToken	= substr($event,$f,1);
				if ($eventToken eq $Break) {
					$retProb *= (1-$retFreq);				# Token specifies fragment is lost (prob = 1-RF)
				} else {
					$retProb *= $retFreq;					# Fragment is retained at retention frequency
					my @contains = split "-", $fragments[$f];
					push @content, @contains;				# Add the markers to the content
				}
			}
			my $stuff = join '', (sort {$a cmp $b} @content);
			$_ = $stuff;
			s/(.)\1+/$1/g;
			$_ .= substr($nullString, 0,2-length($_));
			$subContent{$_}	+= $retProb;
		}
		foreach my $n (keys %subContent) {
			next if ($subContent{$n} <= 0);
			$Content{$k}{$n} = $subContent{$n};
		}
	}
	for my $d (0..$#Distances) {
		my $dist		= $Distances[$d];
	#	print "<FONT COLOR=blue><B>$self->{Name} $dist kb</B></FONT>\n";
		my $breakProb	= 1 - exp(-$lambda * $dist);
	#	print "breakProb = $breakProb ";
		foreach my $k (@pKeys) {
			my $baseProb	= $probability{$k};
	#		print "$k Base prob = $baseProb & ";
			foreach my $i ( @variables ) {						# Now factor in the probabilities for the variable distances
				my $breakToken	= substr($k,1+$i*2,1);
				if ($breakToken eq $Break) {
					$baseProb *= $breakProb;					# Token specifies break
				} else {
					$baseProb *= (1-$breakProb);				# No break, prob = 1 - probabilityOfBreak
				}
	#			print qq(token = "$breakToken", );
			}
	#		print "modified to $baseProb at $dist kb.\n";
			next if ($baseProb == 0);
			foreach my $n ( keys %{$Content{$k}} ) {
				${$outputRef}{$n}[$d] += $Content{$k}{$n} * $baseProb;
			}
		}
	}
#	print "</PRE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub scanMarker {
	my $self		= shift;
	my $markerRef	= shift;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# [0] - width of image
# [1] - name of image file
# [2] - imported global values
sub drawChromosome {
	$hT = (gdTinyFont->height);
	$wT = (gdTinyFont->width);
	my (@mark, @dist);
	my $self		= shift;
	my $width		= shift;
	my $fileName	= shift;
	my $gv			= shift;
	my $dist		= $self->chrParse(\@mark, \@dist);
	my $pad 		= 500;
	my $size		= $dist + $pad * 2;
	my $avgFrag		= ${$gv}{kbPercR} * 100;
	my $bar			= 5;				my $left	= 0;
	my $height		= $bar + $wT*5;	my $right	= 0;
	my $gph			= new GD::Image($width + $left + $right,$height);
	my $white		= $gph->colorAllocate(255,255,255);
	my $red			= $gph->colorAllocate(255,0,0);
	my $green		= $gph->colorAllocate(0,128,0);
	my $blue		= $gph->colorAllocate(0,0,255);
	my $black		= $gph->colorAllocate(0,0,0);
	my $gray		= $gph->colorAllocate(150,150,150);
	my $orange		= $gph->colorAllocate(255,128,0);
	my $col			= $blue;
	$gph->transparent($white);
	$gph->interlaced('true');
	my $fragRight = $left + $width * $avgFrag / $size;
	$gph->filledRectangle($left,0, $fragRight, $bar, $orange);
	$gph->string(gdTinyFont,$fragRight,0,"  Average Fragment Size",$orange);
	
	my $ybar	= $height-$bar;
	my $ytick	= $ybar - 5;
	my $ylabel	= $ytick - 2;
	
	my @sorted = ();
	foreach my $n (keys %{$self->{Markers}}) {
		foreach my $pos (@{$self->{Markers}{$n}}) {
			my $fracPos = $pos / $size;
			push @sorted, "$fracPos\t$n";
		}
	}
	@sorted = sort {$a cmp $b} @sorted;
	my $x = $pad; my $lastx = $width * $x / $size; my $c = $col;
	$gph->filledRectangle($left,$ybar, $lastx, $height, $gray);
	for my $i (0..$#mark) {
		my $x2 = $xPlot = $left + $width * $x / $size;
		$x2 = $lastx + $hT + 1 if ($i > 0 && $x2 < $lastx + $hT);
		$gph->line($x2,$ytick,$xPlot,$ybar,$col);
		$gph->filledRectangle($lastx,$ybar, $xPlot, $height, $c);
		$gph->stringUp(gdTinyFont,$x2 - $hT/2,$ylabel,$mark[$i],$col);
		$lastx = $xPlot;
		if ($dist[$i] eq "V") {
			$c = $green;
			$x += $self->{Token}{VAR};
		} elsif ($dist[$i] eq "U") {
			$c = $white;
			$x += $self->{Token}{UNL};
		} else {
			$x += $dist[$i];
			$c = $col;
		}
	}
	$gph->filledRectangle($lastx,$ybar, $left+$width, $height, $gray);

	open (FILE, ">$fileName") or die "Failure to write to $fileName: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$fileName"><B>$self->{Name}</B><BR>\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This package manages the statistics generated from the virtual chromosomes
package chrStat;
@ISA = (Chromosome);									# May not need it, but make it inherit the Chromosome package
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub new {
	my $class		= shift;
	my $name		= shift;
	my $valRef		= shift;
	my $probRef		= shift;
	my $distRef		= shift;
	my $tokRef		= shift;
	my $self = {};
	@{$self->{Distance}}	= ();
	%{$self->{Probability}}	= ();
	%{$self->{Token}}		= ();
	$self->{Model}			= $name;
	foreach my $n (keys %{$valRef}) {		# Set all the general variables
		$self->{$n} = ${$valRef}{$n};		 
	}
	foreach my $n (keys %{$probRef}) {		# Set the probability arrays
		@{$self->{Probability}{$n}} = ();
		for my $p (0..$#{${$probRef}{$n}}) {
			push @{$self->{Probability}{$n}}, ${$probRef}{$n}[$p];
		}
	}
	@{$self->{Distance}} = ();
	for my $p (0..$#{$distRef}) {
		push @{$self->{Distance}}, ${$distRef}[$p];
	}
	foreach my $n (keys %{$tokRef}) {		# Set the token hash
		$self->{Token}{$n} = ${$tokRef}{$n};
	}
	my $nt = 0; my @An = my @Bn = ();
	my @AB = @{$self->{Probability}{"AB"}};
	my @nn = @{$self->{Probability}{"$nt$nt"}};
	if (exists $self->{Probability}{"A$nt & B$nt"}) {
		@An = @Bn = @{$self->{Probability}{"A$nt & B$nt"}}
	} else {
		@An = @{$self->{Probability}{"A$nt"}};
		@Bn = @{$self->{Probability}{"B$nt"}};
	}
	my $elements = $#{$distRef};
	my $k = "d" . ($elements+1);
	my $num = Probability::addStatistic($name,$elements, pack($k,@{$distRef}), pack($k,@AB), pack($k,@An), pack($k,@Bn), pack($k,@nn));
#	print ">$num ";
	bless $self, $class;
	return $self;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub junk {
	@testx = (1.777,17,45,90,50,55); $t = $#testx+1;
	$p = pack("d" . $t,@testx); $ll = length ($p);
	print qq(Sent "$p" with length $ll. <BR>\n);
	$x = Probability::testThang(5,"Heya", $p, $#testx);
	print "-We got $x<BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This package holds linkage data parsed from an output file
package link;
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub new {
	my $class		= shift;
	my $data		= shift;				# [0] The text line straight from the file
	my $gv			= shift;				# [1] ref to hash of global values
	my $self		= {};
	my $s			= ${$gv}{sep};			# Field separator
	my $LC			= ${$gv}{lodCutoff};
	my @elements	= split "\t", $data;
	$self->{A}		= shift @elements;		# Name of A marker
	$self->{B}		= shift @elements;		# Name of B marker
	$self->{link}	= false;				# Are the markers linked?
	my @first		= split $s, $elements[0];
	foreach my $h (@elements) {				# Look at each linkage hit
		my @d = split $s, $h;
		push @{$self->{prob}}, $d[0];		# -log probability for the hit
		push @{$self->{dPrb}}, $d[1];		# standard error for probability
		push @{$self->{dist}}, $d[2];		# distance in kb
		push @{$self->{dDst}}, $d[3];		# standard error for distance
		push @{$self->{modl}}, $d[4];		# model for this hit
		unless ($d[0] - $first[0] > $LC) {	# Compare each model to the first one (highest prob)
			push @{$self->{nice}}, $#{$self->{modl}};
		}
		
		if ($h !~ /Unl/i) {
			$self->{link} = true;
		}
	}
	bless $self, $class;
	return $self;
}
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
