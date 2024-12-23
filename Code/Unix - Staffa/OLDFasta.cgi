#!/usr/local/bin/perl
BEGIN {
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser. From Steve's Primer3 program
    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
}

print "<html><TITLE>Execution Script</TITLE><BODY BGCOLOR=White>";

$theCommand = "ls -l";
print "<H1>Executing: $theCommand</H1><P><P><PRE>\n";
system "$theCommand";
die;

################################################
# Required Directory Structure and Support Files
# - ParseQuery.pl		Routine for parsing HTML arguments
# - CGI Package			Lincoln Stein's CGI library
# - GD Package			Lincoln's GIF graphics manipulation library
# Data/					Folder at same level as Fasta.cgi
# Data/Alignments/		Folder to contain sequence alignments
# GIF/					Folder at same level as Fasta.cgi for storing temporary GIF files *** CAN BE REDIFINED IN THE DEFINECONSTANTS SUBROUTINE AS $ResLoc
# Data/RepeatLibrary	Fasta Library of genomic repeats

use CGI;
use GD;
require "ParseQuery.pl";

&DEFINECONSTANTS;
&PARSE;

if ($File) {
	print "<html><TITLE>Sequence Alignment</TITLE><BODY BGCOLOR=White>";
	&SHOWALIGNS;
} elsif ($AlignThem ne "") {
	print "<html><TITLE>Internal Library Alignment</TITLE><BODY BGCOLOR=White>";
} else {
	print "<html><TITLE>FASTA Web Interface</TITLE><BODY BGCOLOR=White>";
	&PRINTOPTIONS;
#	print qq(<FORM ACTION=Fasta.cgi><INPUT TYPE=SUBMIT VALUE="Random"><INPUT TYPE=TEXT NAME=MakeRandom SIZE=10> bp</FORM>);
}

if ($othername) {
	$otherstuff = "";
	while (<$othername>) {
		$otherstuff .= $_;
	}
}

undef (@F);
if ($AlignThem ne "") {
	# A special case that will take machine generated files and do FastA alignment on them
	$filename = "Q" . $AlignThem . "$Here.aln"; $libname = "L" . $AlignThem . "$Here.aln";
	open (FILE, "$HOME/Data/Alignments/$filename") or die "Can't read $filename: $!\n";
	while (<FILE>) {
		$F[0] .= $_;
	}
	close FILE;
	open (FILE, "$HOME/Data/Alignments/$libname") or die "Can't read $libname: $!\n";
	while (<FILE>) {
		$F[1] .= $_;
	}
	close FILE;
	&PARSEREF;
} elsif ($filename || $LibCompare) {
	while (<$filename>) {
		$F[0] .= $_;
	}
	if ($libname) {
		while (<$libname>) {
			$F[1] .= $_;
		}
	}
	&PARSEREF;
}


&RANDOMSEQ($MakeRandom) if ($MakeRandom);

# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PARSEGENSCAN {
	($d,$b) = split "Gn.Ex", $otherstuff;
	($a,$d) = split "Predicted", $b;
	@lines = split /[\n\r]/, $a;
	print "<PRE>";
	undef (@Assembled);
	for $i (2..$#lines) {
		next if ($lines[$i] eq "");
		undef (@parts);
		@parts = split " ", $lines[$i];
		($num,$element) = split /\x2E/, $parts[0];
		$element += 0;
		print "<FONT COLOR=blue>>$parts[1]";
		if ($parts[1] =~ "Intr" || $parts[1] =~ "Exon") {
			print " ", $element + 0;
		}
		print " for Transcript $num</FONT>\n";
		$start = $parts[3];
		$start = $parts[4] if ($parts[2] eq "-");
		$count = 0; $Region = substr ($QuerySeq,$start-1, $parts[5]);
		$Region = &REVERSE($Region) if ($parts[2] eq "-");
		while ($count < length($Region)) {
			print substr ($Region, $count, $BlockSize), "\n";
			$count += $BlockSize;
		}
		unless ($parts[1] =~ "PlyA" || $parts[1] =~ "Prom") {
			$Assembled[$num][$element] = $Region;
		}
	}
	for $i (0..$#Assembled) {
		$Region = "";
		for $j (1..$#{$Assembled[$i]}) {
			$Region .= $Assembled[$i][$j];
		}
		next if ($Region eq "");
		print "<FONT COLOR=blue>>Predicted mRNA Transcript $i</FONT>\n";
		$count = 0;
		while ($count < length($Region)) {
			print substr ($Region, $count, $BlockSize), "\n";
			$count += $BlockSize;
		}
	}
	print "</PRE>";
	$Region = "";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Just dumps the contents of a file(s) to a window. Used to show alignments
sub SHOWALIGNS {
	undef (@Files);
	@Files = split "-", $File;
	if ($#Files > 0) {
		print "<FONT COLOR=orange><B>This match represents multiple overlapping alignments.</B></FONT><BR>";
	}
	
	print "<PRE>";
	foreach $i (@Files) {
		open (FILE, "$HOME/Data/Alignments/$i") or die "Can't read $i: $!\n";
		while (<FILE>) {
			print $_;
		}
		close FILE;
		print "<HR>";
	}
	print "</PRE><TABLE>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PRINTOPTIONS {
	print "<TABLE BORDER=1>";
	print "<TR><TD ALIGN=center COLSPAN=3 BGCOLOR=GOLD><B><FONT SIZE=+2>FASTA Web Interface</FONT></B><BR>";
	print "<i><FONT COLOR=olive>Mouse-over any field for help on the functioning of each setting</i>";	
	$colwidth = 250;
	
	print "<TR><TD VALIGN=TOP WIDTH=$colwidth>";
	print $query->start_multipart_form(
								-action=>"Fasta.cgi",
								-method=>POST,
								-name=>"MainForm");
	print "<FONT COLOR=darkgreen SIZE=+1><B>Query file</B></FONT><BR>\n";
	print $query->filefield(	-name=>'uploaded_file',
								-default=>$filename,
								-title=>"A text file containing a single DNA sequence in FastA format. Use the Browse... button to select the file, do not type it in directly.",
								-size=>20),"<BR>\n";
	print "<FONT COLOR=brick SIZE=+1><B>Library file</B></FONT><BR>\n";
	print $query->filefield(	-name=>'uploaded_lib',
								-title=>"A text file containing one or more DNA sequences in FastA format.",
								-default=>$libname,
								-size=>20),"<BR>\n";
	print "<FONT COLOR=orange SIZE=+1><B>Other Files</B></FONT><BR>\n";
	print $query->filefield(	-name=>'other_file',
								-title=>"A text file containing supplemental information. Currently only Genescan output files can be interpreted.",
								-default=>$othername,
								-size=>20),"<BR>\n";
	print $query->submit(		-value=>'Analyze the indicated files',
								-title=>"Run the program with the settings and files indicated above.",
								),"\n";
	print "<TD VALIGN=TOP WIDTH=$colwidth>\n\n";
	print "Chunk size of \n";							
	print $query->textfield(	-name=>'ChunkSize',
								-title=>"This is the maximum length of the query sequence that will be analyzed at a time. Smaller chunks will take longer to analyze, but will allow closely neighboring repeats to be detected.",
								-default=>$ChunkSize,
								-size=>6),"bp,\n";
	print " overlap of \n";							
	print $query->textfield(	-name=>'ChunkOverlap',
								-title=>"Indicates how many bases each chunk should overlap. This value should be small compared to the chunk size.",
								-default=>$ChunkOverlap,
								-size=>3),"bp<BR>\n";
	print $query->checkbox(-name=>'CheckRepeat',
								-title=>"Will run an additional FastA search of your query against a local library of known genomic repeats.",
								-checked=>$CheckRepeat,
								-label=>'',
								-value=>'1'), "\n Check for Repeat Elements<BR>\n";
	print "Find Exons ";
	print $query->checkbox(-name=>'FindExons',
								-checked=>$FindExons,
								-title=>"To use this feature, the query file should contain genomic sequence of the gene of interest, and the Library file should contain the cDNA/EST sequence. Works well only if the genomic sequence and cDNA show a high degree of conservation.",
								-label=>'',
								-value=>'1');
	print " Word Size ";
	print $query->textfield(	-name=>'wordSize',
								-title=>"Genomic and cDNA sequence are compared by searching for 100% matches between 'words', small fragments of DNA. Larger fragments will give greater specificity, smaller fragments will allow more divergent sequences to be compared (but will more likely fail).",
								-default=>$wordSize,
								-size=>2)," bp<BR>\n";
	print "&nbsp;&nbsp;<FONT SIZE=-1> Show Introns ";
	print $query->checkbox(-name=>'ShowIntrons',
								-title=>"Check if you wish a FastA dump of predicted introns as well as exons. This may increase load time if the introns are large.",
								-checked=>$ShowIntrons,
								-label=>'',
								-value=>'1'), " Skip Fasta Search ";
	print $query->checkbox(-name=>'NoFasta',
								-title=>"Normally, the program will perform a fasta search of the predicted exons against the query. Check if you wish to skip this (perhaps because fasta searching is slow at the time).",
								-checked=>$NoFasta,
								-label=>'',
								-value=>'1'), "</FONT><BR>";
	print "Self-compare library ";
	print $query->checkbox(-name=>'LibCompare',
								-title=>"When checked, the program will just compare all library sequences against themselves (no other operations will be performed). A very simple but surprisingly handy algorithm will list the top scoring matches. This feature is a useful compainion to Sequencher, which often does not assemble chimeric clones.",
								-checked=>$LibCompare,
								-label=>'',
								-value=>'1');
	print " Word Size ";
	print $query->textfield(	-name=>'libSize',
								-title=>"The Self-Compare search relies on word matching. Larger words increase the specificity, but may miss matches in poor quality sequence. Smaller words will find matches in poor sequence, but will find more junk as well.",
								-default=>$libSize,
								-size=>2)," bp<BR>\n";
	print "Rank hits by ";
	print $query->popup_menu(	-name=>'SortBy',
								-title=>"The criteria that you want FastA matches to be sorted by. Higher ranked hits will appear first in the text listing. Note that Z-scores are not generated if the library is small, while Smith-Waterman scores are always generated.",
								-values=>[qw/SWscore Zscore PerID/],
								-labels=>{'SWscore'=>'Smith-Waterman', 'Zscore'=>'Score', 'PerID'=>'Percent Identity'},
								-default=>'$SortBy'), "<BR>\n";
	print $query->checkbox(-name=>'ColorType',
								-title=>"Hits will not be color-coded by rank, but will be coded by category - the top 10 library entries will have unique colors assigned to them, all other entries will have a common color.",
								-checked=>$ColorType,
								-label=>'',
								-value=>'1'), "Color Code by Library Entry<BR>";
	print "Cutoff at Category ";
	print $query->popup_menu(	-name=>'Cutoff',
								-title=>"Setting a cutoff will prevent the display of lower-scoring hits. Hits are binned into 10 score categories, shown on the key seen after analysis. A black line will be drawn at the cutoff category when a cutoff is specified.",
								-values=>[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
								-labels=>{'0'=>'None'},
								-default=>'$Cutoff'), "\n";
	print "<TD VALIGN=TOP WIDTH=$colwidth>";
	print $query->checkbox(-name=>'DoCPG',
								-title=>"Perform a CpG dinucleotide search of the Query and display the results graphically.",
								-checked=>$DoCPG,
								-label=>'',
								-value=>'1'), "<B>CpG Analysis</B> &nbsp;";
	print $query->textfield(	-name=>'minBinSize',
								-title=>"The search will tally the number of CpGs found in a continuously scanning window. Smaller windows provide more precision, but will also generate higher jitter in the graphs, whereas larger windows provide smoother graphs and more of a broad picture. RESULTS ARE HIGHLY DEPENDANT ON WINDOW SIZE, EXPERIMENT WITH VALUES.",
								-default=>$minBinSize,
								-size=>4)," bp window<BR>\n";
	print $query->checkbox(-name=>'DoPalindrome',
								-title=>"Check the Query sequence for palindromes and hairpins, and display the results graphically. This search is realtively computationally intensive.",
								-checked=>$DoPalindrome,
								-label=>'',
								-value=>'1'), "Palindrome Analysis<BR>";
	
	print "&nbsp;&nbsp;<FONT SIZE=-1> Look for Palindromes of ";
	print $query->textfield(	-name=>'MinPalSize',
								-title=>"Minimum size palindromes to graph. Do not recomend dropping below 4!",
								-default=>$MinPalSize,
								-size=>2),"-\n";
	print $query->textfield(	-name=>'MaxPalSize',
								-title=>"Maximum size palindromes to graph. Don't be shy, keep it large - it will not adversley affect computational time.",
								-default=>$MaxPalSize,
								-size=>2)," bp<BR>\n";
	print "&nbsp;&nbsp; Minimum Hairpin size ";
	print $query->textfield(	-name=>'MinHairpin',
								-title=>"Minimum size of one hairpin arm to consider.",
								-default=>$MinHairpin,
								-size=>2)," with<BR>\n";
	print "&nbsp;&nbsp; loop length of ";
	print $query->textfield(	-name=>'MinLoop',
								-title=>"Minimum size of hairpin loop (the sequence between the palindromic arms).",
								-default=>$MinLoop,
								-size=>2),"-\n";
	print $query->textfield(	-name=>'MaxLoop',
								-title=>"Maximum size of hairpin loop. Do not make this too large, it will greatly increase computation time.",
								-default=>$MaxLoop,
								-size=>2)," bp<BR>\n";
	print "&nbsp;&nbsp; Analyze 4-mers ";
	print $query->checkbox(-name=>'FourFrequency',
								-title=>"Generate a table of relative frequency of each (16 total) 4-mer. I think this is not working currently...",
								-checked=>$FourFrequency,
								-label=>'',
								-value=>'1'), "</FONT><BR>";
	print "Image Width ";
	print $query->textfield(	-name=>'width',
								-title=>"The size of the main graphic in pixels. Larger sizes give you more spatial resolution, but may be a pain to view with your browser.",
								-default=>$width,
								-size=>4)," pixels<BR>\n";
	print $query->endform;
	print "</TABLE>";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub RANDOMSEQ {
	srand ( time() ^ ($$ + ($$<<15)) );
	if ($TestRandom) {
		$Sequences[0] = "";
		for $i (0..$_[0]-1) {
			$Sequences[0] .= $RandNucs[rand(4)];
		}
	} else {
		print "<PRE>";
		print ">Random ",$_[0]," bp";
		for $i (0..$_[0]-1) {
			print "\n" unless ($i % $BlockSize);
			print $RandNucs[rand(4)];
		}
		print "</PRE>";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADFASTA {
	undef (@a); undef (@b); undef (@junk);
	@junk = split /[\x00\x01\x05]/, $_[0];					# Desperately trying to eliminate crap brought in with the resource fork by internet explorer.
	$j = -1; $firstChar = "DUMMY";
	if ($#junk == 0) {
		@a = split /[\r\n]/, $junk[0];
	} else {
		while ($rc{$firstChar} eq "" && $j <= $#junk) {				# Once a line is found with a >, does the next line start with an allowed char?
			$j++;
			while ($j <= $#junk && substr($junk[$j],0,1) ne ">") {	# Find first line of junk that starts with >
				$j++;
			}
			@a = split /[\r\n]/, $junk[$j];					# @a hopefully contains the data fork
			$firstChar = uc( substr ($a[1],0,1));			# the first character of the second line - should be a nucleotide!
		}
	}
#	&LOADERROROUT;										# De-comment to dump de-bugging info - used when trying to figure out how to elimintate resource fork.
	return if ($j > $#junk || $#a < 1);					# did not find a ">", not a FASTA file
	@b = split ">", $a[0]; $a[0] = ">" . $b[$#b];		# Try to eliminate leading junk from mac files. This is needed for internet explorer
	$last = $#a;
	while ( $a[$last] =~ /[^$Allowed]/i || $a[$last] eq "") {		# find the last line that contains only allowed characters
		$last--;
	}
	@a = @a[0..$last];							# this eliminates trailing junk 
	return @a;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Screen dump variables for LoadFasta
sub LOADERROROUT {
	print "First line is $a[0]<BR>\n";
	print "First Char at $j is $firstChar<BR>\n";
	print "Last Line is $a[$#a]<BR>\n";
	print "Size of a is $#a. J is $j, size of junk is $#junk<BR>\n";
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
sub PARSEFASTA {
	$Tag = $_[0]; $Seq = ""; $num = 0;
	undef (@Names); undef (@Sequences);
	for $i (1..$#_) {
		next unless ($_[$i]);
		$n = substr( $_[$i] , 0 ,1);
		$Seq .= $_[$i] unless ($n eq ">"); # attach next line of sequence unless it has a FASTA tag
		if ($i == $#_ or $n eq ">") {
			while (substr($Tag, -1, 1) eq " ") { chop $Tag; }	# Remove trailing spaces
			push @Names, $Tag;
			push @Sequences, uc($Seq);
			$num++;
			$Tag = $_[$i];
			$Seq = "";
		}
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PALINDROME {
	$_ = $QueryName; s/\>//; $QueryName = $_;	# Eliminate ">" in names
	$end = $NumberOfBP - 2;	# second to last character
	undef (@PalHits); undef (@SimpleHits); undef (@Hairpins); undef (%Hit);

	######################################
	# Scan along all bases in the sequence
	$PalTime = time;
	for ($Pos = 0; $Pos <= ($end - $MinReach); $Pos++) {
		
		######################
		# Simple Repeat Search
		if (($Pos + $MinSimple) < $end) {
			$k = 0;
			while ($k < $SimpSize) {			# Scan for simple repeats (may be confused as palindromes eg ATATAT)
				$k += $SimpStep;
				$numRepeats = 1;
				$basicRept = substr ($QuerySeq, $Pos, $k);
				while ($basicRept eq (substr ($QuerySeq, $Pos + ($k*$numRepeats), $k)) ) {
					$numRepeats++;
				}
				last if ($numRepeats >= $SimpLimit[$k]);	# Stop if the number of repeat units is at least that specified by SimpLimit
			}
			if ($numRepeats >= $SimpLimit[$k]) {			# A simple repeat has been found
				push @{$SimpleHits[$k*$numRepeats][$k]}, $Pos;	# Add to @SimpleHits
#print "$RepeatNames[$k] Repeat at <B>$Pos</B> $numRepeats x $basicRept = ", substr ($QuerySeq, $Pos, $k*$numRepeats) , "<BR>\n";
				$Pos += $k * $numRepeats - 1;				# Move counter one ahead
				next;
			}
		}
		###############################
		# Palindrome & Hairpin Search
		for $loop (0,$MinLoop..$MaxLoop) {
			next if ($Hit{$Pos - $loop});		# The region has already been identified as harboring a palindrome
			$reach = $MaxReach;
			$reach = $Pos if (($Pos - $reach) < 0);
			$reach = $end - $Pos if (($reach + $Pos + 1) > $end);
			last if ($reach < $MinReach);
			$res = -1;
			for $k (0..$reach) {
				$a = substr ($QuerySeq, $Pos - $k - $loop, 1);	# The left base
				$b = substr ($QuerySeq, $Pos + $k + 1, 1);		# The mirror imaged right base
				last if ($rc{$a} ne $b);						# stop searching if they are not complements
				$res++;
			}
			
			if ($loop == 0) {
				if ($res >= $MinReach) {
					$siz = ($res+1)*2;
					push @{$PalHits[$siz]}, $Pos;
					$excludedDistance = $res - $MinHairpin;
					$excludedDistance = 0 if ($excludedDistance < 0);
					for ($k = $Pos - $excludedDistance; $k <= $Pos; $k++) { $Hit{$k} = $siz - ($Pos - $k); }	# 
					last;
				}
			} else {
				if ($res >= $MinHairpin) {
					push @{$Hairpins[($res+1)*2 + $loop][$loop]}, $Pos - $res - $loop;
#print "Hairpin at <B>$Pos</B> with loop $loop = ", &COLORSEQ(substr ($QuerySeq, $Pos - $res-$loop, ($res+1)*2 + $loop),$loop) , "<BR>\n";
					$excludedDistance = $res - $MinHairpin;
					$excludedDistance = 0 if ($excludedDistance < 0);
					for ($k = $Pos - $excludedDistance - $loop; $k <= ($Pos - $loop); $k++) { $Hit{$k} = $siz - ($Pos - $k); }	# 
					last;
				}
			}
		}
	}

	$max = $maxPal = $#PalHits;
	$maxHair = ($#Hairpins);
	$SimpleSpace = 0;		# Number of bases occupied by simple repeats 
	for $k (0..$#SimpleHits) {
		for $rptSiz (0..$#{$SimpleHits[$k]}) {
			for $j (0..$#{$SimpleHits[$k][$rptSiz]}) {
				$SimpleSpace += $k;
			}
		}
	}
	
	$max = $#SimpleHits if ($#SimpleHits > $max);
	$max = int($maxHair/2) if ($maxHair/2 > $max);
	
	$Pheight = 200; $Ptop = 10; $Pbottom = 5;
	$Pbottom = $wT * 7 if ($ShowPalTicks);
	
	$gph = new GD::Image($width + $left + $right,$Pheight + $Ptop + $Pbottom);
	$white = $gph->colorAllocate(255,255,255);
	$red = $gph->colorAllocate(255,0,0);
	$green = $gph->colorAllocate(0,128,0);
	$blue = $gph->colorAllocate(0,0,255);
	$ltblue = $gph->colorAllocate(100,100,255);
	$black = $gph->colorAllocate(0,0,0);
	$orange = $gph->colorAllocate(255,128,0);
	$ltgreen = $gph->colorAllocate(0,255,0);
	$darkred = $gph->colorAllocate(128,0,0);
	$gray50 = $gph->colorAllocate(128,128,128);
	$gray75 = $gph->colorAllocate(192,192,192);
	$gph->transparent($white);
	$gph->interlaced('true');
	
	@Cols = ($red, $green, $ltblue);
	
	$nummark = $width / ($hT + 2);
	$vrat = $Pheight / $max;
	
	$gph->line($left,$top,$left,$Pheight+$Ptop,$black);
	$gph->line($width+$left,$Pheight+$Ptop,$left,$Pheight+$Ptop,$black); $k = 0;
	
	$m = "Size of feature in bp"; $HalfHeight = ($wT * length($m) / 2) + ($Pheight + $Ptop) / 2;
	$gph->stringUp(gdTinyFont,0,$HalfHeight,$m,$gray50);
	for ($k = 0; $k <= $max; $k ++) {	# Draw horiz grid lines
		$y  = $Pheight+$Ptop - $k * $vrat; $kc++;
		$col = $black;
		$gph->line($left,$y,$left-3,$y,$col);
		$gph->string(gdTinyFont,$left-4-($wT*length($k)),$y-($hT/2),$k,$black) unless ($k % 5);
		next if ($k % 2);
		if ($k % 10) {
			$gph->dashedLine($left,$y,$width+$left,$y,$col);
		} else {
			$gph->line($left,$y,$width+$left,$y,$col);
		}
	}
	
	$CircleColor = $gray75;  $CircleDia = 10;
	$dx = int (0.5 + $wT/2); $dy = int (0.5 + $hT/2);
	##################################
	# Draw Repeat Icons
	$RepeatColor = $orange;
	for ($k = $max; $k >= $MinPalSize; $k--) {
		$y = $Pheight+$Ptop - $k * $vrat;
		for $rptSiz (0..$#{$SimpleHits[$k]}) {
			for $j (0..$#{$SimpleHits[$k][$rptSiz]}) {
				$x = $SimpleHits[$k][$rptSiz][$j] * $ratio + $left;
				for ($fill=2; $fill <=$CircleDia; $fill +=2) {				# Other fill routines are unsatisfactory
					$gph->arc($x,$y,$fill,$fill,0,360,$CircleColor);
				}
				$gph->string(gdTinyFont,$x - $dx,$y - $dy,$rptSiz,$blue);
			}
		}
	}
	
	##################################
	# Draw Palindrome Marks
	for ($k = $max; $k >= $MinPalSize; $k--) {
		if ($#{$PalHits[$k]} >= 0) {
			$y = $Pheight+$Ptop - $k * $vrat;
			$col = $Cols[ ($k/2) % ($#Cols + 1) ];
			for $j (0..$#{$PalHits[$k]}) {
				$x = $PalHits[$k][$j] * $ratio + $left;
				$gph->line($x,$y,$x,$Pheight+$Ptop,$col);
			}
		}
	}
	
	$CircleColor = $orange;
	#########################
	# Draw Hairpin Icons
	for ($k = $#Hairpins; $k >= $MinPalSize; $k--) {
		if ($#{$Hairpins[$k]} >= 0) {
			for $lp (0..$#{$Hairpins[$k]}) {
				if ($#{$Hairpins[$k][$lp]} >= 0) {
					$y = $Pheight+$Ptop - (($k-$lp)/2) * $vrat;
					for $j (0..$#{$Hairpins[$k][$lp]}) {
						$x = $Hairpins[$k][$lp][$j] * $ratio + $left;
						$gph->filledRectangle($x-$dx-1,$y-$dy-1,$x+$dx+1,$y+$dy+1,$CircleColor);
						$gph->string(gdTinyFont,$x - $dx,$y - $dy,"H",$black);
					}
				}
			}
		}
	}
	
	
	##################################
	# Add x-axis ticks and values:
	$maxTicks = $width / ($hT +2);	# The maximum number of x-axis values that can legibily be displayed
	$tickStep = int ($NumberOfBP / $maxTicks);				# The spacing that such ticks would have
	@NiceTicks = (1,2,5,10);
	$factor = 1; $nt = 0;
	while ($NiceTicks[$#NiceTicks] * $factor < $tickStep) { $factor *= 10; }	# Find the order of magnitude that encompases the tickstep
	while ($NiceTicks[$nt] * $factor < $tickStep) { $nt++;}						# Find the appropriate "nicetick" that suits tickstep

	for ($k = 0; $k <= $NumberOfBP; $k += $NiceTicks[$nt] * $factor) {
		$x = $k * $ratio + $left;
		$gph->line($x,$Pheight+$Ptop,$x,$Pheight+$Ptop + 3,$black);
		$gph->stringUp(gdTinyFont,$x-($hT/2),$Pheight+$Ptop + 6 + $wT*length($k),$k,$black);
	}
	&FOURPLOT if($FourFrequency);
	
	############################
	# Save file and print tables
	$gphname = "Palgraph$$" . $Here . ".tmp";
	open (FILE, ">$ResLoc/$gphname") or die "sub PALINDROME: Failure to write to $gphname: $!\n";
	print FILE $gph->gif;
	close FILE;
	
	$PalTime = time - $PalTime;
	print qq(<img src="$ResWeb/$gphname"><BR>\n);
}
# - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + 
# Plot local frequency of 4-base palindromes
sub FOURPLOT {
	$bintime = time;
	$basesPerPixel = $NumberOfBP / $width;			# number of bases each pixel represents
	$pixelReach = 1;										# number of pixels on either side of current pixel to consider
	$minBinSize = $BinSize;
	while ( ($pixelReach*2 + 1) * $basesPerPixel < $minBinSize) { $pixelReach++; }
	$fourBinSize = ($pixelReach*2 + 1) * $basesPerPixel;

	undef (@Bin); undef (%FourHit);
	for ($palsize = 4; $palsize <= $#PalHits; $palsize +=2) {
		next if ($#{$PalHits[$palsize]} < 0);
		for $k (0..$#{$PalHits[$palsize]}) {
			$target = int($PalHits[$palsize][$k] / $basesPerPixel);	# The bin/pixel the current hit falls in
			$low = $target - $pixelReach;
			$low = 0 if ($low < 0);
			$pn = substr ($QuerySeq, $PalHits[$palsize][$k]-1, 4);
			$FourHit{ $pn }++;	# Counts the specific kinds of palindromes
			$FourTot{ $pn }++;	# Keeps track of total hits
			for $j ($low..($target+$pixelReach)) {
				$Bin[$j]++;
			}
		}
	}
	
	$binMax = 0; $CGbinMin = 100000000;
	for $k ($pixelReach..($width-$pixelReach)) {
		$binMax = $Bin[$k] if ($binMax < $Bin[$k]);
		$binMin = $Bin[$k] if ($binMin > $Bin[$k]);
	}
	$Divisor = 16; $NiceStep = 0.05; $NiceMultiple = 1 / $NiceStep;
	$BinExpected = $fourBinSize / $Divisor;													# The number of 4mers we expect in a bin
	$topRatio = $NiceStep + int ($NiceMultiple * $binMax / $BinExpected) / $NiceMultiple;	# An upper observed/expected ratio that encompases the highest value
	$bottomRatio = int ($NiceMultiple * $binMin / $BinExpected) / $NiceMultiple;			# As above, but the lower value
	$binBaseline = $bottomRatio * $BinExpected;												# The lowest number of hits to be plotted
	$binRat = (4 * $vrat) / ($BinExpected * ($topRatio - $bottomRatio));					# vertical plotting ratio
	
	$Draw{$bottomRatio} = "Y"; $Draw{$topRatio} = "Y";
	$NumNiceTicks = $NiceStep;
	while ( (4 * $vrat) / (($topRatio - $bottomRatio) / $NumNiceTicks) < ($hT + 1)) {	# Make sure there are not too many ticks for this scale
		$NumNiceTicks += $NiceStep; 
	}
	
	$y  = $Pheight+$Ptop - int(0.5 + ($BinExpected - $binBaseline) * $binRat);
	$gph->line($left,$y,$width+$left,$y,$gray50);
	for ($k = $bottomRatio; $k <= $topRatio; $k += $NumNiceTicks) {
		$y = $Pheight+$Ptop - int(0.5 + $binRat * ($k * $BinExpected - $binBaseline));
		$gph->line($width+$left,$y,$width+$left+3,$y,$black);
		$gph->string(gdTinyFont,$width+$left+5,$y-($hT/2),$k,$black);
	}
	
	$yOld 	= $Pheight+$Ptop - int(0.5 + ($Bin[$pixelReach] - $binBaseline) * $binRat);
	$xOld = $pixelReach * $basesPerPixel * $ratio + $left;
	for $k ($pixelReach..($width-$pixelReach)) {
		$x = $k * $basesPerPixel * $ratio + $left;
		$y =	$Pheight+$Ptop - int(0.5 + ($Bin[$k] - $binBaseline) * $binRat);
		$gph->line($x,$y,$xOld,$yOld,$black);
		$xOld = $x; $yOld = $y;
	}
	$bintime = time - $bintime;
}
# - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + 
sub PALSUMMARY {
	$BinSize = int (0.5 + $BinSize);	# Clean up for printing
	$BinExpected = int (0.5 + $BinExpected);
	print "\n\n<TABLE BORDER=1><TR>";
	print "<TD><B>Size</TD><TD><B>Observed</TD><TD><B>Expect</TD><TD><B>Obs/Exp</TD></TR>\n";
	$Summary .= "\n<TR><TD>$QueryName</TD><TD>" . $NumberOfBP . "</TD>\n";
	
	for ($k = $MinPalSize; $k <= $MaxPalSize; $k += 2) {
		$expected = $NumberOfBP - $SimpleSpace;		# Don't count simple repeat bases when calculating expected frequencies
		$found = 0;
		$res = $k/2 -1;
		for $j (0..$res) { $expected /= 4; }
		for $j ($k..$#PalHits) { $found += $#{$PalHits[$j]} + 1; }	# Each larger palindrome also represents an example of the current one
		$rat = $found / $expected;
		if ($rat > 100) {
			$rat = int(0.5 + $rat);
		} else {
			$rat = int (0.5 + $rat *100) / 100;
		}
		$output = sprintf("<TD>%d</TD><TD>%d <FONT COLOR=GRAY>(%d)</TD><TD>%d</TD><TD>%s</TD>" , $k , $found, $#{$PalHits[$k]} + 1 , int ($expected + 0.5) , $rat);
		$Summary .= "<TD>$rat</TD>";
		$TotalFound[$k] += $found;
		next if ($found < 1);
		print "<TR>$output</TR>\n";
	}
	$TotalBases += $NumberOfBP;
	$TotalSimple += $SimpleSpace;
	print "<TR><TD COLSPAN=4><FONT SIZE=-3>";
	print "P[4] bin size = $BinSize bp (expect $BinExpected per bin)<BR>\n" if ($FourFrequency);
	print "Total time = $PalTime sec<BR>\n";
	print "$SimpleSpace bases occupied by simple repeats.</FONT>\n";
	print "</TABLE>\n\n";
}
# - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + 
#	[0...15] = array of values		[n-1] = size of bar	[n] = total number of bases
sub FOURDISTRIBUTION {
	@theKeys = keys %FourHit;
	$max = 0; $Bar = ""; $numicons = 30; undef (@ToSort);
	$expect = $NumberOfBP / 256; $prettyExpect = int (0.5 + $expect);
	for $k (0..$numicons) { $Bar .= "#"; }
	
	foreach $k (@theKeys) {
		$max = $FourHit{$k} if $FourHit{$k} > $max;
		push @ToSort, sprintf ("%10.d\t$k",  $FourHit{$k});
	}
	@ToSort = sort {$b cmp $a} @ToSort;
	print "<B><FONT COLOR=orange SIZE=+1>4-mer summary. Expect $prettyExpect for each palindrome.</FONT></B><BR>\n";
	print "<PRE>";
	foreach $k (@ToSort) {
		($a, $b) = split "\t", $k;
		printf ("%s <B>%-${numicons}.${numicons}s</B>  %3.0d    <FONT COLOR=gray>(%3.2f)</FONT>\n", 
					$b, substr($Bar,0, int (0.5 + $FourHit{$b} * $numicons / $max)) , $FourHit{$b}, $FourHit{$b} / $expect);
	}
	print "</PRE>";
}
# - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + 
# [0] = index for $Sequences
sub BESTHITS {
	print "<TABLE BORDER=1><TR><TH>Position</TH><TH>Size</TH><TH>Sequence</TH></TR>";
	if ($maxHair) {
		$numTop = 20;
		print "<TR><TD COLSPAN=3 ALIGN=CENTER BGCOLOR=khaki><FONT COLOR=forestgreen SIZE=+2><B>Hairpins over length $numTop in $QueryName</FONT></TD></TR>\n";
		for ($k = $maxHair; $k >= $numTop; $k--) {
			for $loop ($MinLoop..$MaxLoop) {
				if ($#{$Hairpins[$k][$loop]} >= 0) {
					for $j (0..$#{$Hairpins[$k][$loop]}) {
						print "<TR><TD><B>$Hairpins[$k][$loop][$j]</TD><TD>$k</TD><TD><PRE><B>";
						print &COLORSEQ(substr ($QuerySeq, $Hairpins[$k][$loop][$j], $k),$loop), "</TD></TR>";
					}
				}
			}
		}
	}
	$numTop = 15;
	print "<TR><TD COLSPAN=3 ALIGN=CENTER BGCOLOR=khaki><FONT COLOR=NAVY SIZE=+2><B>Palindromes over length $numTop in $QueryName</FONT></TD></TR>\n";
	for ($k = $max; $k >= $numTop; $k--) {
		if ($#{$PalHits[$k]} >= 0) {
			for $j (0..$#{$PalHits[$k]}) {
				print "<TR><TD><B>$PalHits[$k][$j]</TD><TD>$k</TD><TD><PRE><B>";
				print &COLORSEQ(substr ($QuerySeq, $PalHits[$k][$j] - ($k/2 - 1), $k)), "</TD></TR>";
			}
		}
	}
	$numTop = 25;
	print "<TR><TD COLSPAN=3 ALIGN=CENTER BGCOLOR=khaki><FONT COLOR=BRICK SIZE=+2><B>Simple Repeats over length $numTop in $QueryName</FONT></TD></TR>\n";
	for ($k = $max; $k >= $numTop; $k--) {
		for $rptSiz (0..$#{$SimpleHits[$k]}) {
			for $j (0..$#{$SimpleHits[$k][$rptSiz]}) {
				print "<TR><TD><B>$SimpleHits[$k][$rptSiz][$j]</TD><TD>$k</TD><TD><PRE><B>";
				print &COLORSEQ(substr ($QuerySeq, $SimpleHits[$k][$rptSiz][$j], $k)), "</TD></TR>";
			}
		}
	}
	print "</TABLE>\n";
}
# - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + 
#	[0] = the sequence to be analyzed	[1] = the number of bases in the loop of a hairpin
sub COLORSEQ {
	$a = $_[0];
	$l = (length($a) - $_[1]) / 2;	# number of bases on either side of loop
	$out = ""; $c = 0;
	while ($a ne "") {
		$_ = chop $a; $c++;
		if ($_[1] && ($c > $l) && ($c <= $l + $_[1])) {
			$out = "<FONT COLOR=purple SIZE=-2>$_</FONT>" . $out;
		}elsif (/a/i) {
			$out = "<FONT COLOR=green>A</FONT>" . $out;
		} elsif (/c/i) {	
			$out = "<FONT COLOR=blue>C</FONT>" . $out;
		} elsif (/t/i) {	
			$out = "<FONT COLOR=red>T</FONT>" . $out;
		} else {	
			$out = $_ . $out;
		}
		
	}
	return $out;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Extracts data values from hit records previously stored in dataline
sub PARSEDATALINE {
	foreach $dd (@dataTypes) {		# Clear variables
		$$dd = "";
	}
	undef (@data); undef (@Filenames);
	@data = split "\t", $_[0];
	foreach $dd (@data) {
		($variable, $value) = split "=", $dd;
#		print "$variable = $value<BR>\n";
		$$variable = $value;
	}
#	@Filenames = split "-", $Files;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# Reassembles the data line and returns it
sub SETDATALINE {
	$theData = "";
	foreach $dd (@dataTypes) {
		$theData .= $dd . "=" . $$dd . "\t";
	}
#	$Files = join "-", @Filenames;
#	$theData .= $Files;
	return $theData;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub SETCOLOR {
	$white		= $gph->colorAllocate(255,255,255);
	$black		= $gph->colorAllocate(0,0,0);
	$gph->transparent($white);
	$gph->interlaced('true');
	$gray50		= $gph->colorAllocate(128,128,128);
	$gray75		= $gph->colorAllocate(192,192,192);
	$fog		= $gph->colorAllocate(92,186,201);
	$watercress	= $gph->colorAllocate(187,212,26);
	$fern		= $gph->colorAllocate(92,163,31);
	$marigold	= $gph->colorAllocate(255,204,23);
	$orange		= $gph->colorAllocate(255,145,23);
	$melon		= $gph->colorAllocate(255,117,23);
	$salmon		= $gph->colorAllocate(234,97,87);
	$carnation	= $gph->colorAllocate(217,69,107);
	$apple		= $gph->colorAllocate(234,31,28);
	$orchid		= $gph->colorAllocate(186,87,194);
	$tan		= $gph->colorAllocate(191,133,74);
	$dirt		= $gph->colorAllocate(140,97,54);
	$lemmon		= $gph->colorAllocate(255,235,23);
	
	undef (@spect); # This is the color spectrum used by TNG RH
	$spect[0]	= $white;
	$spect[1]	= $gph->colorAllocate(214,214,214);
	$spect[2]	= $gph->colorAllocate(189,189,189);
	$spect[3]	= $gph->colorAllocate(156,156,156);
	$spect[4]	= $gph->colorAllocate(0,255,206);
	$spect[5]	= $gph->colorAllocate(0,255,90);
	$spect[6]	= $gph->colorAllocate(33,255,0);
	$spect[7]	= $gph->colorAllocate(148,255,0);
	$spect[8]	= $gph->colorAllocate(255,239,0);
	$spect[9]	= $gph->colorAllocate(255,123,0);
	$spect[10]	= $gph->colorAllocate(255,0,0);

	
	$red		= $gph->colorAllocate(255,0,0);
	$green		= $gph->colorAllocate(0,128,0);
	$yellow		= $gph->colorAllocate(255,255,0);
	$blue		= $gph->colorAllocate(0,0,255);
	$lightblue		= $gph->colorAllocate(128,128,255);
	for ($k=20; $k <= 100; $k += 5) {
		$j = int(255 * (100-$k) / 100);
		$blue{$k} = $gph->colorAllocate($j,$j,255);
	}
	$blue		= $gph->colorAllocate(0,0,255);
	$darkred	= $gph->colorAllocate(128,0,0);
	
###	@score =	(0,			30,			50,			100,		200,		300,			350,		400,		450,		500);
#	@ScoreCol	= ($white,	$gray75,	$gray50,	$fog,		$fern,		$watercress,	$marigold,	$melon,		$apple,		$orchid,	$green);
	@ScoreCol	= @spect;
	@RepeatCol	= ($white,	$blue{20},	$blue{30},	$blue{40},	$blue{50},	$blue{60},		$blue{70},	$blue{80},	$blue{90},	$blue{95},	$blue{100});
	if ($ColorType) {
		for $i (0..$#ScoreCol) { $ScoreCol[$i] = $gray50; }
		@TagCol = ($marigold, $lightblue, $melon, $red, $lemmon, $orchid, $fern, $watercress, $blue, $tan, $dirt);
		$TagCol[$TagCounter] = $dirt;
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub WRITELIB {
	undef (%Lookup); undef (@Reverse);
	open (FILE, ">$HOME/Data/Library$Here") or die "sub WRITELIB: Failure to write to $HOME/Data/Library$Here: $!\n";
	for $i (0..$#Names) {
		$Names[$i] = substr ($Names[$i], 0, 20) unless ($GoExon > 1);
		print FILE "$Names[$i]\n";
		print FILE "$Sequences[$i]\n";
		$Lookup{"$Names[$i]"} = $i;
		$Lookup{"$Names[$i] RC"} = $i;
		print FILE "$Names[$i] RC\n";
		$Reverse[$i] = &REVERSE($Sequences[$i]);
		print FILE "$Reverse[$i]\n";
	}
	close FILE;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub REVERSE {
	$d = $_[0]; $rev = "";
	while (length($d) > 0) {
		$n = chop $d;
		if ($rc{$n}) {
			$rev .= $rc{$n};
		} else {
			$rev .= $n;
		}
	}
	return $rev;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub DOFASTA {
	if ($FindExons ne "") {
		&WORDCHECK;
		$SortBy = "PerID";
	}
	
	if ($NumberOfBP <= $ChunkSize) {
		$ChunkSize = $NumberOfBP;												# Single Chunk Analysis
		$NumberOfChunks = 1;
	} else {
		$NumberOfChunks = int($NumberOfBP / $ChunkSize) + 1;
		$x = int ( ($NumberOfBP + $ChunkOverlap *($NumberOfChunks - 1) ) / $NumberOfChunks) + 1;
		while ($x > $ChunkSize) {			# The above calculations should work in almost all cases. However, there will be exceptions, especially if ChunkSize is close to Overlap
			$NumberOfChunks++; $Outofhand++;
			$x = int ( ($NumberOfBP + $ChunkOverlap *($NumberOfChunks - 1) ) / $NumberOfChunks) + 1;
			if ($Outofhand > 100) {
				print "<B>Uh... I'm having problems calculating the appropriate number of chunks to break up the query into.<BR>";
				print "Perhaps you've chosen too small a chunk size, or too large an overlap. The overlap should be small in comparison to the chunksize.<BR>";
				print "(for example, chunksize = 5000 and overlap = 40). Try again with different values...<BR>";
				die;
			}
		}
		$ChunkSize = $x;
	}
	
#	elsif ( ($ChunkSize * 2 - $ChunkOverlap) > $NumberOfBP) {
#		$ChunkSize = int (($NumberOfBP - $ChunkOverlap) / 2) + $ChunkOverlap;	# Two Chunk Analysis
#		$NumberOfChunks = 2;
#	} else {
#		$NumberOfChunks = ($NumberOfBP - 2 * ($ChunkSize - $ChunkOverlap/2));	# Remove the chunks from the edge (which have only one overlap each)
#		$NumberOfChunks = $NumberOfChunks / ($ChunkSize - $ChunkOverlap);		# Find the number of chunks in the middle of the sequence
#		$NumberOfChunks = int($NumberOfChunks) + 1 + 2;							# Round up (int + 1) and then add in the two edge chunks
#		$ChunkSize = int ( ($NumberOfBP + $ChunkOverlap *($NumberOfChunks - 1) ) / $NumberOfChunks);
#	}
	
	@rank = @{$Rank{$SortBy}};
	
	print "<!"; system "rm $HOME/Data/Alignments/*$Here.aln"; print ">";
	$fastaTime = $processTime = 0;

	undef (@Alignments); 
	$Start = 0;
	
	while ($Start + $ChunkOverlap < $NumberOfBP) {
		undef (@output);
		$chunkcount++;
		open (FILE, ">$HOME/Data/Query$Here") or die "Failure to write to Query: $!\n";
		printf FILE ("%s from base %d to %d\n", $QueryName, $Start+1, $Start+$ChunkSize);
		# NEED TO MAKE SURE THAT THE LAST BASE GETS READ!! - May have to add a fudge to ChunkSize
		print FILE substr($QuerySeq, $Start, $ChunkSize);
		print FILE "\n";
		close FILE;
		
		$ti = time; 
		# -h		Turn off Histogram
		# -z 0		Turn off Statistics
		# -Q		Quiet mode - also allows use of b/d modifiers
		# -b 60		Number of simularity Scores to be shown
		# -d 60		Number of best alignments to be shown
		# -n		Forces query to be treated as DNA
		$NumHits = 60;
		$output[0] = ` /usr/local/etc/httpd/cgi-bin/pagelab/fasta3  -A -n -Q -H -z 0 -b $NumHits -d $NumHits $HOME/Data/Query$Here $HOME/Data/Library$Here` if ($libname);
		$output[1] = ` /usr/local/etc/httpd/cgi-bin/pagelab/fasta3  -A -n -Q -H -z 0 -b $NumHits -d $NumHits $HOME/Data/Query$Here $HOME/Data/RepeatLibrary` if ($CheckRepeat);
		$fastaTime += (time - $ti);
		
#		print "<PRE>$output[0]</PRE>";
		
		$ti = time;
		for $k (0..$#output) {
			next if ($output[$k] eq "");
			undef (@parse); undef (@temp);
			@parse = split ">>", $output[$k];				# split output into alignments
			if ($parse[0] !~ "Scan time") {
				print "<TABLE><TR><TD BGCOLOR=yellow><FONT COLOR=RED><B>Probable FastA error while analyzing chunk $chunkcount</B></FONT> <i>Text \"Scan time\" not found in output. Full output follows:</i></TABLE>\n";
				print "<pre><FONT COLOR=navy>$output[$k]</FONT></pre>";
			}
			@temp = split "\n\n\n", $parse[$#parse];		# Get rid of trailing data (statistics and such)
			$parse[$#parse] = $temp[0];
			for $parseCount (1..$#parse) {
				undef (@content); $DataLine=""; undef (@stuff); undef (@temp);
				@content = split "\n", $parse[$parseCount];	# split alignments into lines
				@temp = split "  ", $content[0];
				################################################
				# Store the following variables in the data line:
				# initn, init1, opt, Z-score, expect
				$Name = $temp[0];
				$DataLine = "Name=" . $Name;
				$_ = $content[1]; s/\://g; s/\(\)//g; s/Z\-/Z/; s/rev-comp//;
				@stuff = split " ", $_;							# split up the score line into var/values
				for ($k = 0; $k <= $#stuff; $k += 2) {
					$DataLine .= "\t" . $stuff[$k] . "=" . $stuff[$k+1];
				}
				$content[2] =~ /score:\s*(\d+);\s*([\d\.]+)% identity in (\d+) nt overlap/;
				$DataLine .= "\tSWscore=" . $1;
				$DataLine .= "\tPerID=" . $2;
				$DataLine .= "\tOverlap=" . $3;
				$filename = "Alignment" . $alignCounter++ . "$Here.aln";	# The filename that will be used for the alignment files
				$firstmatch = index($content[6], ":");	 				# Find the first occurance of a match in the first line
				$MatchStart	= &FINDPOSITION($firstmatch,	$content[4],			$content[5])			+ $Start;
				$LibStart	= &FINDPOSITION($firstmatch,	$content[8],			$content[7]);
				$lastLine = $#content;
				while ($content[$lastLine] !~ ":") { $lastLine--; }		# Find the last line with a match
				$lastmatch = rindex($content[$lastLine], ":");	 		# Find the last occurance of a match in the last line
				$MatchEnd	= &FINDPOSITION($lastmatch,		$content[$lastLine-2],	$content[$lastLine-1])	+ $Start;
				$AllowSecond = "Y";
				$LibEnd		= &FINDPOSITION($lastmatch,		$content[$lastLine+2],	$content[$lastLine+1]);
				if ($LibEnd < $LibStart) {								# Occasionally the lower match will not have number tags on it
					$AllowSecond = "";
					$firstlastmatch = index($content[$lastLine], ":");	# Find the **First** occurance of a match in the last line
					$remainder = $lastmatch - $firstlastmatch + 1;		# This is the number of matched bases on the last line
					$upline = length($content[$lastLine-5]) -1;			# Should correspond to last base of the second-to-last DNA match line.
					$LibEnd		= &FINDPOSITION($upline,		$content[$lastLine-4],	$content[$lastLine-5]);
					$LibEnd += $remainder;
#					print "$upline, $lastmatch - $firstlastmatch = $remainder<BR>";
				}
				
				$DataLine .= "\tMatchStart=$MatchStart\tMatchEnd=$MatchEnd";
				$DataLine .= "\tLibStart=$LibStart\tLibEnd=$LibEnd\tFiles=$filename" . "-";
				push @Alignments, $DataLine;			# Store the alignment info
# print "------<BR>\n"; $_ = $DataLine; s/\t/\<BR\>/g; print "$_<BR>\n";
				
				##########################
				# Write the alignment to a file. Include FastA formated sequence for the aligned stretches.
				open (FILE, ">$HOME/Data/Alignments/$filename") or die "sub DOFASTA: Can't write to $filename: $!\n";
				print FILE "<FONT SIZE=+2 COLOR=NAVY><B>$Name - Add $Start to all values.</B></FONT>\n";
				print FILE $parse[$parseCount];
				$_ = $QueryName; s/\>//; $n = $_;
				print FILE "\n>Matched Region for $n\n";
				$count = 0; $Region = substr ($QuerySeq,$MatchStart-1, $MatchEnd-$MatchStart+1);
				while ($count < length($Region)) {
					print FILE substr ($Region, $count, $BlockSize), "\n";
					$count += $BlockSize;
				}
				
				$index = $Lookup{">$Name"};						# Find the index for the match
				if ($index eq "") {
					print FILE "\n<FONT COLOR=orange><i>FastA format for library alignment not available</i></FONT>\n";
				} else {
					$_ = $Names[$index]; s/\>//; $n = $_;
					if (($LibEnd - $LibStart + 1) == length($Sequences[$index]) ) {			# The whole region was matched
						$Region = $Sequences[$index];
					} else {
						$a = $LibStart - 1;
						if ($Name =~ "RC") {
							$a = length($Sequences[$index]) - $LibEnd;
						}
						$Region = substr ($Sequences[$index],$a, $LibEnd-$LibStart+1);
					}
					$count = 0; 
					print FILE "\n>Matched Region for $n\n";
					while ($count < length($Region)) {
						print FILE substr ($Region, $count, $BlockSize), "\n";
						$count += $BlockSize;
					}
					print FILE "\n<B><FONT COLOR=green>The whole subject was matched.</B></FONT>\n" if (($LibEnd - $LibStart + 1) == length($Sequences[$index]) );
					print FILE "<FONT COLOR=brick>The normal complement is shown for the match.</FONT>\n" if ($Name =~ "RC");
				}
				close FILE;
			}
		}
		$Start += ($ChunkSize - $ChunkOverlap);
		$processTime += (time - $ti);
	}
	
	###############################
	# Check to see if any alignments of the same thing overlap (possible with chunked analysis).
	undef (@ToSort); undef (@ScoreSort); undef (@Final); undef (@temp);
	for $i (0..($#Alignments-1)) {
		&PARSEDATALINE($Alignments[$i]);
		$ThisName	= $Name;				# Variables to hold needed info from the queried alignment
		$ThisStart	= $MatchStart;
		$ThisEnd	= $MatchEnd;
		$ThisRank	= $$SortBy;
		$TheseFiles	= $Files;
		$concatenate = "";
		for $j (($i+1)..$#Alignments) {
			&PARSEDATALINE($Alignments[$j]);
			if ($ThisName eq $Name) {												# Two matches with the same name
				if ( $ThisStart <= $MatchEnd && $ThisStart >= $MatchStart ) {		# The start of a is within b
					if ($ThisEnd > $MatchEnd) {										# The match of a extends further than b
						$MatchEnd = $ThisEnd;
					}
					$concatenate = "Y";
				}
				if ( $ThisEnd <= $MatchEnd && $ThisEnd >= $MatchStart ) {			# The end of a is within b
					if ($ThisStart < $MatchStart) {									# The match of a extends further than b
						$MatchStart = $ThisStart;
					}
					$concatenate = "Y";
				}
			}
			if ($concatenate) {														# Overlapping match was found
#print "Concatenating $Name $MatchStart - $MatchEnd<BR>";
				$$SortBy = $ThisRank if ($$SortBy < $ThisRank);						# Set the rank to the higher ranked hit
				$Files = $TheseFiles . $Files;										# Add the alignment matches from a
				$Alignments[$j] = &SETDATALINE;										# Reset the dataline in question
				last;
			}
		}
		push @temp, $Alignments[$i] unless ($concatenate);
	}
	push @temp, $Alignments[$#Alignments];
	
	###############################
	# Eliminate alignments with scores below cutoff
	for $i (0..$#temp) {
		&PARSEDATALINE($temp[$i]);
		if ($$SortBy >= $rank[$Cutoff]) {			# Check to see that the cutoff was ok
			push @Final, $temp[$i];
			push @ToSort, sprintf ("%10.d\t%d",  $MatchStart, $#Final);
			push @ScoreSort, sprintf ("%10.d\t%d",  $$SortBy, $#Final);
		} else {
			push @CutOffed, $temp[$i];
#			print "<FONT COLOR=blue SIZE=-1>$temp[$i]</FONT><BR>\n";
		}
	}
	@ToSort = sort {$a cmp $b} @ToSort;
	@ScoreSort = sort {$b cmp $a} @ScoreSort;
	
	
	###############################
	# Find just the best matches for the exons
	if ($FindExons ne "") {	
		undef (@ScoreSort);
		for $k (0..$#ExonNames) {
			$n = $ExonNames[$k];
			$maxSW=0; $loc = -1;
			for $j (0..$#Final) {
				&PARSEDATALINE($Final[$j]);
				next unless (">$Name" eq $n || ">$Name" eq "$n RC");
				if ($SWscore > $maxSW) {
					$loc = $j;
					$maxSW = $SWscore;
				}
			}
			if ($loc > -1) {
				push @ScoreSort, "null\t$loc";
			} else {
				print "<FONT COLOR=red>No FastA matches found for $ExonNames[$k].";
				if ($ExonNames[$k] =~ "Exon") {
					if (length($Exons[$k]) < 50) {
						print " This could be due to its small size.";
					} else {
						print " <B>This is disturbing - a previously identified exon should have been located.</B>";
					}
				}
				print "</FONT><BR>\n";
			}
		}
		@ToSort = @ScoreSort;
	}
	
	if ($ColorType) {
		$TagCounter = 0;
		$MaxTag = 10;		# Maximum number of classes 
		for $i (0..$#ScoreSort) {
			($dummy, $k) = split "\t", $ScoreSort[$i];
			&PARSEDATALINE($Final[$k]);
			@stripped = split " RC", $Name;			# Get rid of RC suffixes
			if ($ColorTag{$stripped[0]} eq "") {	# A tag has not been assigned to this hit (name) yet
				$ColorTag{$stripped[0]} = $TagCounter;
				if ($TagCounter < $MaxTag) {	# The first $MaxTag classes will be assigned unique tags, the rest will have a generic one, with index of $MaxTag
					$ColorName[$TagCounter] = $stripped[0];
					$TagCounter++;
				}
			}
		}
		$ColorName[$TagCounter] = "All Others";
	}

}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub PARSEREF {
	undef (@temp); undef (@Loaded); $rec = {};
	%rc = (	"A" => "T", C => G, T => A, G => C, );
	@types = ("Query" , "Library");
	for $filenum (0..1) {
		next if ($F[$filenum] eq "");
		@Loaded = &LOADFASTA($F[$filenum]);
		if ($#Loaded < 0) {
			$t = $filename;
			$t = $libname if ($filenum == 1);
			if ($t) {
						$parsed = join "\n", @Loaded;
print <<EOF;
<P><FONT COLOR=RED SIZE=+1><B>$t does not seem to be a FastA file.</B></FONT><BR>
<P>This annoyance generally occurs because some bizzare characters got included during a copy-and-paste operation. 
Certain programs written by Big Brother Bill will inapropriately load the file's resource fork as well as the data fork. 
This has the result of confusing the parsing algorithm. 
The best solution so far seems to be to <FONT COLOR=RED><B>use Netscape</B></FONT>. 
If you must still use Explorer, then try:<BR><BR>
<OL>
<LI>Copy the sequence only (not the name tag) out of your original FastA file.<BR>
<LI> Paste it into EditSeq.<BR>
<LI> Select all (command-A) and copy it.<BR>
<LI> Go to BBEdit, open a new document, and paste the sequence in.<BR>
<LI> Add a new FastA name tag (e.g. >MySequence)to the top, and save the file.<BR>
<LI> Try again with the new file.<BR>
</OL>
<P><B>Here is how I cleaned up the file:</B><FONT COLOR=navy>
<PRE>$parsed</PRE>
</FONT><P><B>This is the contents of the raw file as loaded by the browser:</B><FONT COLOR=navy>
<PRE>$F[$filenum]</PRE>
</FONT><FONT COLOR=white>
EOF
				die;
			}
			next;
		}
		
		$numLines = $#Loaded + 1;
		&PARSEFASTA(@Loaded);
		$num = $#Names + 1; $numbp = 0;
		for $k (0..$#Names) {
			$numbp += length ($Sequences[$k]);
		}
		if ($filenum == 0) {	# Store Query data in data structures - may need to save it as chunks
			$QueryName = $Names[0];
			$QuerySeq  = $Sequences[0];
			$NumberOfBP = length ($QuerySeq);
			$out = &TIDYBROWSER($filename);
			print "<FONT COLOR=brown><B>Query file <FONT COLOR=blue>$out</FONT> contains $numLines lines and $NumberOfBP bp.</B></FONT><BR>\n";
			if ($num > 1) {
				print "<FONT COLOR=RED><B>WARNING:</B></FONT> The query file contained $num entries. Only the first entry has been used: ", $numbp - $NumberOfBP, "bp were ignored!<BR>\n";
			}
		} else {
			&WRITELIB;
			$out = &TIDYBROWSER($libname);
			print "<FONT COLOR=brown><B>Library file <FONT COLOR=blue>$out</FONT> contains $numLines lines and $num entries, total sequence $numbp bp.</B></FONT><BR>\n";
		}
		undef (@Loaded);
	}
	undef (@F);
	
	if ($LibCompare) {
		&COMPARELIBRARY;
		die;
	}
	
	$ratio = $width / $NumberOfBP;
	$left = 30;	$right = 30;
	
	if (($libname ne "" && $NoFasta eq "") || $CheckRepeat) {
		&DOFASTA;
		&TOPSTUFF;
	} else {
		$_ = $QueryName; s/>//; 
		print "<CENTER><FONT SIZE=+1 COLOR=green><B>Results for <FONT COLOR=brick>$_</FONT> Analysis</B></FONT></CENTER><BR>\n";
	}
	
	#########################################
	# Calculate Aesthetic X-axis tick spacing
	$maxTicks = $width / ($hT + 5);	# The maximum number of x-axis values that can legibily be displayed
	$tickStep = int ($NumberOfBP / $maxTicks);				# The spacing that such ticks would have
	@NiceTicks = (1,2,5,10);
	$factor = 1; $nt = 0;
	while ($NiceTicks[$#NiceTicks] * $factor < $tickStep) { $factor *= 10; }	# Find the order of magnitude that encompases the tickstep
	while ($NiceTicks[$nt] * $factor < $tickStep) { $nt++;}						# Find the appropriate "nicetick" that suits tickstep

	&PALINDROME if ($DoPalindrome);
	&CPGSEARCH if ($DoCPG);
	
	$FeatureHeight = $top = 1;			# If no FastA analysis done, then don't leave any extra graphical space
	$height = $FeatureHeight + ($wT * 10);
	if ($libname || $CheckRepeat) {			
		$FeatureHeight = 9;				# Maximum feature height
		$top = 10;
		$dnaAT = 20;					# Y position of DNA info
		$height += $dnaAT + ($hT * 10); # Space to print the alignment names
	}
	
	if ($ColorType) {
		$w = 120;
		$h = $hT * ($#TagCol + 1);
		$gph = new GD::Image($w,$h);
		&SETCOLOR;
		$y = 0;
		for $k (0..$#TagCol) {
			$m = "$SortBy $rank[$k-1] - $rank[$k]";
			$gph->string(gdTinyFont,0,$y,$ColorName[$k],$TagCol[$k]);
			$y += $hT;
		}
		$gphname = "Taggraph$$" . $Here . ".tmp";
		open (FILE, ">$ResLoc/$gphname") or die "Failure to write to $gphname: $!\n";
		print FILE $gph->gif;
		close FILE;
		print qq(<TD VALIGN=top><img src="$ResWeb/$gphname">\n);
	}
	
	print qq(</TABLE><A NAME="TOP"></A>\n);
	print qq(\n<MAP NAME="align">\n);
	
	$gph = new GD::Image($width + $left + $right,$height);	
	&SETCOLOR;
	
	
	##################################
	# Add x-axis ticks and values:
	$yOffset = $dnaAT + $FeatureHeight + $top + 1;
	for ($k = 0; $k <= $NumberOfBP; $k += $NiceTicks[$nt] * $factor) {
		$x = $k * $ratio + $left;
		$gph->line($x,$yOffset,$x,$yOffset + 3,$black);
		$gph->stringUp(gdTinyFont,$x-($hT/2),$yOffset + 6 + $wT*length($k),$k,$black);
	}
	$gph->filledRectangle($left,$dnaAT+$top-1,$width+$left,$dnaAT+$top+1,$black);	# Draw bar to represent sequence

	###############################################################
	# Draw graphical hits from worst to best
	for ($sorting=$#ToSort; $sorting >= 0; $sorting--) {
		($dummy, $k) = split "\t", $ToSort[$sorting];
		&PARSEDATALINE($Final[$k]);
		$x1 = int (0.5 + $MatchStart * $ratio + $left);
		$x2 = int (0.5 + $MatchEnd * $ratio + $left);
		if ($Name =~ " RC") { 		#Reverse Complement
			$y1 = $dnaAT + $top+2;
			$y2 = $dnaAT + $FeatureHeight + $top;
		} else {
			$y1 = $dnaAT - $FeatureHeight + $top;
			$y2 = $dnaAT + $top-2;
		}
		$col = 0;
		for $l (1..$#rank) {
			if ( $$SortBy < $rank[$l]) {
				$col = $l;
				last;
			}
		}
		$col = $#rank + 1 if ($$SortBy >= $rank[$#rank]);
		
		$_ = $Name; s/\'//g; $out = $_; $textcol = $black;
		if ($Name =~ "Repeat") {
			undef (@temp);
			$y1++; $y2--;
			@temp = split " ", $out;
			$out = $temp[1];
		} else {
			$col = $ScoreCol[$col];
			if ($ColorType) {
				@stripped = split " RC", $Name;
				$col = $TagCol[ $ColorTag{$stripped[0]} ];
			}
		}
		
		print qq(<AREA SHAPE=rect COORDS="$x1,$y1,$x2,$y2" );
		print qq(ONMOUSEOVER='document.OUT.tx.value="$out : );
		print "Score=$Zscore, " if ($Zscore);
		print qq(SWScore=$SWscore, %ID=$PerID"' );
		print qq(ONMOUSEOUT='document.OUT.tx.value=""' );
		print qq(HREF="Fasta.cgi?File=$Files" target=_blank>\n);
		$gph->filledRectangle($x1,$y1,$x2,$y2,$col);
	}
	
	###############################################################
	# Print names of hits below the DNA line (best to worst)
	for $sorting (0..$#ToSort) {
		($dummy, $k) = split "\t", $ToSort[$sorting];
		&PARSEDATALINE($Final[$k]);
		$_ = $Name; s/\'//g; $out = $_; $textcol = $black;
		$x1 = int (0.5 + $MatchStart * $ratio + $left);
		$x2 = int (0.5 + $MatchEnd * $ratio + $left);
		$y = $yOffset + $wT * 8;
		if ($Name =~ "Repeat") {
			undef (@temp);
			@temp = split " ", $out;
			$out = $temp[1];
			$textcol = $blue;
		}
		$x3 = $x1 + $wT * length($out);		# The right coordinate of the imagemap box
		$x3 = $x2 if ($x2 > $x3);			# If the drawn rectangle is larger than the text, set the image map borders to it
		do {
			$tracking = 0;
			for $l (($x1-1)..($x3+1)) {
				$index = $gph->getPixel($l,$y + $hT/2);
				$index2 = $gph->getPixel($l,$y+1);
				unless ($index == $white && $index2 == $white) {
					$tracking = 1;
					$y += $hT + 1;
					last;
				}
			}
		} while ($tracking);
		$y3 = $y + $hT;
		print qq(<AREA SHAPE=rect COORDS="$x1,$y,$x3,$y3" );
		print qq(ONMOUSEOVER='document.OUT.tx.value="$out : );
		print "Score=$Zscore, " if ($Zscore);
		print qq(SWScore=$SWscore, %ID=$PerID"' );
		print qq(ONMOUSEOUT='document.OUT.tx.value=""' );
		print qq(HREF="Fasta.cgi?File=$Files" target=_blank>\n);
		$gph->filledRectangle($x1,$y,$x2,$y + $hT,$yellow);
		$gph->string(gdTinyFont,$x1,$y,$out,$textcol);
	}
	print "</MAP>\n";

	$gphname = "Hitgraph$$" . $Here . ".tmp";
	open (FILE, ">$ResLoc/$gphname") or die "Failure to write to $gphname: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$ResWeb/$gphname" ismap usemap="#align" border=0>\n);
	
	if ($FindExons) {
		print "<PRE>";
		$cnt = 0;
		for $k (0..$#Exons) {
			if ($ExonCount[$cnt] == $k) {
				print "</PRE><FONT SIZE+=4 COLOR=GREEN><B>Predicted Exon Structure for $mRNAnames[$cnt]:</B></B></FONT><PRE>";
				$cnt++;
			}
			print "<FONT COLOR=BLUE><B>$ExonNames[$k] $ExonMod[$k]</B></FONT>\n";
			$count = 0; 
			while ($count < length($Exons[$k])) {
				print substr ($Exons[$k], $count, $BlockSize), "\n";
				$count += $BlockSize;
			}
		}
		$cnt = 0;
		for $k (0..$#Introns) {
			if ($ExonCount[$cnt] == $k) {
#				print "</PRE><FONT SIZE+=4 COLOR=GREEN><B>Predicted Intron Structure for $mRNAnames[$cnt]:</B></FONT><PRE>";
				$cnt++;
			}
			next if ($IntronNames[$k] eq "");
			print "<FONT COLOR=BLUE><B>$IntronNames[$k]</B></FONT>";
			print " <FONT COLOR=red>5' site not GT</FONT>" if (substr ($Introns[$k], 0, 2) !~ m/GT/i);
			print " <FONT COLOR=red>3' site not AG</FONT>" if (substr ($Introns[$k], -2, 2) !~ m/AG/i);
			print "\n";
			$count = 0; 
			while ($count < length($Introns[$k])) {
				print substr ($Introns[$k], $count, $BlockSize), "\n";
				$count += $BlockSize;
			}
		}
		print "</PRE>";
	}
	&PARSEGENSCAN if ($othername);
	
	&FASTAHITS if ($libname);
	
	print "<TABLE><TR>\n";
	if ($DoPalindrome) {
		print "<TD VALIGN=TOP>\n";
		&BESTHITS;
		print "<TD VALIGN=TOP>\n";
		&PALSUMMARY;
		if ($FourFrequency) {
			print "<TD VALIGN=TOP>\n";
			&FOURDISTRIBUTION;
		}

	}
	print "</TABLE>\n";
	

	
	print "<i>FastA analysis time was $fastaTime seconds. $processTime seconds were spent analyzing the output.</i><BR>\n" if ($libname);
	print "<i>$wordTime</i><BR>\n" if ($FindExons);
	if ($DoCPG) {
		print "<BR><i>CpG analysis time was $CGtime seconds.</i><BR>\n";
		$AverageRatio = int (0.5 + $AverageRatio * 100) / 100;
		printf ("<FONT COLOR=orange>Average CG/GC Ratio was %2.2f.</FONT><BR>\n", $AverageRatio);
	}
	print "<i>$citation</i><BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub TOPSTUFF {
	$what = "chunk";
	$what .= "s" if ($NumberOfChunks > 1);
	$_ = $QueryName; s/\>//;
	print "<TABLE><TR><TD ALIGN=CENTER>";
	printf ("<CENTER><FONT SIZE=+2 COLOR=navy><B>Analyzing <FONT COLOR=blue>%s</FONT> in %d %s of %d bp",
						$_, $NumberOfChunks, $what, $ChunkSize);
	print ", $ChunkOverlap bp overlap" if ($NumberOfChunks > 1);
	print ".</B></FONT></CENTER><BR>\n";
	print qq(<BR><FORM NAME="OUT">\n);
	print "<FONT SIZE=+1 COLOR=darkorange><B>Move mouse over a feature for scoring info, click feature for alignment</B></FONT><BR>\n";
	print qq(<CENTER><input name=tx size=80></CENTER>\n);
	print qq (</FORM>\n);
	
	$w = 120;
	$h = $hT * ($#rank + 1);
	$gph = new GD::Image($w,$h);
	&SETCOLOR;
	$y = 0;
	for $k (1..$#rank) {
		$m = "$SortBy $rank[$k-1] - $rank[$k]";
		$gph->string(gdTinyFont,0,$y,$m,$ScoreCol[$k]);
#		$gph->string(gdTinyFont,$w-$wT*3,$y,"Rpt",$RepeatCol[$k]) if ($CheckRepeat);
		$y += $hT;
		if ($k == $Cutoff) {
			$gph->line(0,$y-1,$w,$y-1,$black);
		}
	}
	$gph->string(gdTinyFont,0,$y,"$SortBy $rank[$#rank] +",$ScoreCol[$#rank+1]);
	$gph->string(gdTinyFont,$w-$wT*3,$y,"Rpt",$RepeatCol[$#rank+1]) if ($CheckRepeat);


	$gphname = "Keygraph$$" . $Here . ".tmp";
	open (FILE, ">$ResLoc/$gphname") or die "Failure to write to $gphname: $!\n";
	print FILE $gph->gif;
	close FILE;
	system "chmod 666 $ResLoc/$gphname";
	print qq(<TD VALIGN=top><img src="$ResWeb/$gphname">\n);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub FASTAHITS {
	if ($#ScoreSort >= 0) {		# At least one hit was found
		print "<PRE><FONT COLOR=brick>";
		printf (qq(%20.20s %7.7s  %7.7s  %5.5s  %8.8s  %6.6s\n),
			"Sequence", "Score", "SWScore", "%ID", "Position", "Length");
		print "</FONT>";
		for $sorting (0..$#ScoreSort) {
			($dummy, $k) = split "\t", $ScoreSort[$sorting];
			&PARSEDATALINE($Final[$k]);
			printf (qq(<A HREF="Fasta.cgi?File=$Files" target=_blank>%20.20s</A> %7.0f  %7.0f  %5.1f  %8.8s  %6.6s\n),
				$Name, $Zscore>0 ? $Zscore : "", $SWscore, $PerID, $MatchStart, $MatchEnd-$MatchStart+1);
			$AllFiles .= $Files;
		}
		
		print "</PRE>\n";

		if ($AlignThem ne "") {
			$File = $AllFiles;
			print "<HR>\n";
			&SHOWALIGNS;
		}
	} else {
		print "<BR><FONT COLOR=red SIZE=+1><B>No hits were found for this sequence</B></FONT><BR>\n";
		print "<B>Last Fasta Parse:</B><PRE><FONT COLOR=brick>";
		print "$output[0]";
		print "</PRE>\n";
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# This routine will determine the nucleotide number of a given base in the FASTA output
# [0] = The string index of the base in question
# [1] = The string containing the numerical position tags above the bases
# [2] = The string contain the bases themselves
sub FINDPOSITION {
	undef (@temp);
	$gaps = 0; $tracking = $_[0];
	if (substr($_[1], $tracking, 1) eq " ") {							# The position above the base does not include part of a numerical tag
		while ($tracking > 0 && substr($_[1], $tracking, 1) eq " ") {	# Scan left for a number
			$tracking--;
			$gaps-- if (substr($_[2], $tracking, 1) eq "-");			# Keep track of gaps in sequence
		}
	}
	if ($tracking == 0) {												# Looked all the way left and did not find a number value
		$tracking = $_[0]; $gaps = 0;									# Return to original base position
		while (substr($_[1], $tracking, 1) eq " ") {					# Scan right
			$tracking++;
			$gaps++ if (substr($_[2], $tracking, 1) eq "-");			# Keep track of gaps in sequence
		}
	}
	while (substr($_[1], $tracking+1, 1) =~ /\d/) {						# Now find the start (last character) of the number
		$tracking++;
		$gaps++ if (substr($_[2], $tracking, 1) eq "-");
	}
	@temp = split " ", substr($_[1], $tracking-7);						# Break up line into the numbers
	$c = 0;
#	while ($c < $#temp && $temp[$c] == 0) { $c++; }
	if ($AllowSecond eq "" && ($tracking > 90 or $tracking < 0 or $temp[$c] <= 0)) {
		print qq(<FONT COLOR=red>Position ID failure for <A HREF="Fasta.cgi?File=$filename" TARGET=_blank>$Name</A>: </FONT>\n);
		print "Tracking = $tracking, Gaps = $gaps, nearest ID is \"$temp[$c]\". Parsed Lines are:<BR>\n";
		print "<PRE><FONT COLOR=NAVY>$_[1]\n$_[2]</FONT><\PRE>";
	}
	return ($temp[$c] + $_[0] - $tracking + $gaps);						# Calculate the exact nucleotide position of the base in question
}
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
sub CPGSEARCH {
	$CGtime = time;
	
	$halfBin = int($minBinSize/2);
	$CG = $GC = 0;
	for $Pos (0..($halfBin*2)) {			# Initialize counters for sequence prior to main loop
		$rn = substr ($QuerySeq, $Pos, 2);
		$CG++ if ($rn eq "CG");
		$GC++ if ($rn eq "GC");
	}
	
	$GlobalCG = $LocalCG = $GlobalCGRat = $LocalCGRat = 0; $AverageRatio = 0;
	$LocalCGMin = $LocalCGRatMin = 1000000;
	undef (@CG); undef (@CGRat);
	$xMin = int(($halfBin+1)*$ratio);
	$xMax = int(($NumberOfBP-$halfBin)*$ratio);
	for ($Pos = $halfBin+1; $Pos <= $NumberOfBP-$halfBin; $Pos++) {	# Main loop
		$ln = substr ($QuerySeq, $Pos-$halfBin-1, 2);	# Decrement when a pair leaves the window
		$CG-- if ($ln eq "CG");
		$GC-- if ($ln eq "GC");
		$rn = substr ($QuerySeq, $Pos+$halfBin, 2);	# Increment when a pair comes into the window
		$CG++ if ($rn eq "CG");
		$GC++ if ($rn eq "GC");
		if ($GC == 0) {
			$CGRat = 0;
		} else {
			$CGRat = $CG / $GC;
		}
		$AverageRatio += $CGRat;
		
		$GlobalCGRat =	$CGRat if ($CGRat > $GlobalCGRat);
		$LocalCGRat =	$CGRat if ($CGRat > $LocalCGRat);
		$LocalCGRatMin = $CGRat if ($CGRat < $LocalCGRatMin);
		$GlobalCG =		$CG if ($CG > $GlobalCG);
		$LocalCG =		$CG if ($CG > $LocalCG);
		$LocalCGMin =	$CG if ($CG < $LocalCGMin);
		if (($x = int($Pos*$ratio)) < int(($Pos+1)*$ratio) ) {	# We are about to move into a new pixel, dump data and move on
			$CG[$x] = $LocalCG;				# We take the highest value for the data 'under' this pixel
			$CGRat[$x] = $LocalCGRat;
			$CGMin[$x] = $LocalCGMin;
			$CGRatMin[$x] = $LocalCGRatMin;
			$LocalCGRat = $LocalCG = 0;
			$LocalCGMin = $LocalCGRatMin = 1000000;
		}
	}
	unless (($x = int($Pos*$ratio)) < int(($Pos+1)*$ratio) ) {	# In case it was not added above, which is likely
		$CG[$x] = $LocalCG;				# We take the highest value for the data 'under' this pixel
		$CGRat[$x] = $LocalCGRat;
		$CGMin[$x] = $LocalCGMin;
		$CGRatMin[$x] = $LocalCGRatMin;
	}
	$AverageRatio /= ($NumberOfBP-2*$halfBin);
	
	$Ctop = $Cbottom = 10; $Cheight = 60;
	$gph = new GD::Image($width + $left + $right,$Cheight*2 + $Ctop*2 + $Cbottom);
	$white = $gph->colorAllocate(255,255,255);
	$blue = $gph->colorAllocate(0,0,255);
	$green = $gph->colorAllocate(0,255,0);
	$orange = $gph->colorAllocate(255,128,0);
	$black = $gph->colorAllocate(0,0,0);
	$ltred = $gph->colorAllocate(255,128,128);
	$gray75 = $gph->colorAllocate(192,192,192);
	$gph->transparent($white);
	$gph->interlaced('true');
			
#	$gph->filledRectangle(0,$Cheight + 1.5*$Ctop,$width+$left+$right,$Cheight + 1.5*$Ctop+1,$gray75);
	#####################
	# Plot CG/GC ratio Y axis ticks
	$Ratcol = $orange;
	$m = "CpG/GpC"; $HalfHeight = ($wT * length($m) / 2) + ($Cheight + $Cbottom + $Ctop) / 2;
	$gph->stringUp(gdTinyFont,0,$HalfHeight,$m,$Ratcol);
	$step = 0.1; $NumNiceTicks = 1; $k = 0;
	$RatRat = $Cheight / $GlobalCGRat; 		# CpG vertical plotting ratio
	while ( $NumNiceTicks * $step * $RatRat  < $hT + 2) { $NumNiceTicks++; }
	while ($k <= $GlobalCGRat) {
		$y =	$Cheight + $Ctop - int(0.5 + $k * $RatRat);
		$gph->line($left,$y,$width+$left,$y,$gray75);
		$gph->line($left-3,$y,$left,$y,$black);
		$d = " " . $k;
		$gph->string(gdTinyFont,$left-5 - $wT * length($d),$y-($hT/2),$d,$Ratcol);
		$k += $step * $NumNiceTicks;
	}
	for ($k = 0; $k <= $NumberOfBP; $k += $NiceTicks[$nt] * $factor) {	# Add X-grid
		$x = $k * $ratio + $left;
		$gph->line($x,$y,$x,$Cheight + $Ctop,$gray75);
	}

	
	
	#####################
	# Plot Absolute CpG Y axis ticks
	$CGcol = $blue;
	$m = "CpG per window"; $HalfHeight = ($wT * length($m) / 2) + ($Cheight + $Cbottom + $Ctop) / 2 + $Cheight + $Ctop;
	$gph->stringUp(gdTinyFont,$width + $left + $right - $hT,$HalfHeight,$m,$CGcol);
	$CGbinRat = $Cheight / $GlobalCG; 	# CpG vertical plotting ratio
	$scale = 100; $NumNiceTicks = 1; $k = 0;
	while ( $NumNiceTicks * $CGbinRat < $hT + 2) { $NumNiceTicks++; }
	while ($k <= $GlobalCG) {	
		$y =	$Cheight*2 + $Ctop*2 - int(0.5 + $k * $CGbinRat);
		$gph->line($width+$left,$y,$width+$left+3,$y,$black);
		$gph->line($left,$y,$width+$left,$y,$gray75);
		$gph->string(gdTinyFont,$width+$left+5,$y-($hT/2),$k,$CGcol);
		$k += $NumNiceTicks;
	}
	for ($k = 0; $k <= $NumberOfBP; $k += $NiceTicks[$nt] * $factor) {	# Add X-grid
		$x = $k * $ratio + $left;
		$gph->line($x,$y,$x,($Cheight + $Ctop)*2,$gray75);
	}
	

	$tag = "Window Size: ";
	$gph->string(gdTinyFont,$left,$Cheight + $Ctop+1,$tag,$green);
	$gph->filledRectangle($left + $wT*length($tag),$Cheight + $Ctop+2, $left + ($wT*length($tag)) + ($halfBin * 2 + 1) * $ratio, $Cheight + $Ctop + $hT, $green);
	$ymax = $ratmax = 99999;
	$yOld =		$Cheight + $Ctop - int(0.5 + $CG[$xMin] * $CGbinRat);
	$ratOld =	$Cheight + $Ctop - int(0.5 + $CGRat[$xMin] * $RatRat);
	$xOld =		$xMin + $left;
	for $xp ($xMin..$xMax) {
		$x = $xp + $left;
		$y =	$Cheight + $Ctop - int(0.5 + $CG[$xp] * $CGbinRat) + $Cheight + $Ctop;
		$y2 =	$Cheight + $Ctop - int(0.5 + $CGMin[$xp] * $CGbinRat) + $Cheight + $Ctop;
		$rat =	$Cheight + $Ctop - int(0.5 + $CGRat[$xp] * $RatRat);
		$rat2 =	$Cheight + $Ctop - int(0.5 + $CGRatMin[$xp] * $RatRat);
		$gph->line($x,$y2,$x,$y,$CGcol);
		$gph->line($x,$rat2,$x,$rat,$Ratcol);
#		$gph->line($x,$y,$xOld,$yOld,$CGcol);
#		$gph->line($x,$rat,$xOld,$ratOld,$Ratcol);
#		$xOld = $x; $yOld = $y; $ratOld = $rat;
	}
	
	$gphname = "CpGgraph$$" . $Here . ".tmp";
	open (FILE, ">$ResLoc/$gphname") or die "Failure to write to $gphname: $!\n";
	print FILE $gph->gif;
	close FILE;
	print qq(<img src="$ResWeb/$gphname"><BR>\n);
	$CGtime = time - $CGtime;
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
# This appears to be a bit too memory intensive
sub MAKEWORDHASH {
	$ti = time;
	undef (%forwardHash); undef (%reverseHash);
	for $k (0..($NumberOfBP - $wordSize)) {
		push @{$forwardHash{ substr($QuerySeq, $k, $wordSize) }}, $k;
	}
	for $k (0..(length($QueryRev) - $wordSize)) {
		push @{$reverseHash{ substr($QueryRev, $k, $wordSize) }}, $k;
	}
	$ti = time - $ti;
	$wordTime .= "$ti sec to generate hashtables,  ";
	print "Hashtables created in $ti seconds.<BR>\n";
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub HASHSEARCH {
	$ti = time; $checking = 0;
	$theHash = \%forwardHash; $genomic = \$QuerySeq;
	while ($checking < 2) {
		undef (@AlignStart); undef (@AlignEnd); undef (@QueryStart); undef (@QueryEnd);
		for $k (0..(length($Sequences[$j]) - $wordSize)) {
			undef (@a);
			@a =  @{$$theHash{ substr($Sequences[$j], $k, $wordSize) }};
			if ($#a == 0) {				# Only consider unique matches
				push @AlignStart,	$k;
				push @AlignEnd,		$k    + $wordSize-1;
				push @QueryStart,	$a[0];
				push @QueryEnd,		$a[0] + $wordSize-1;
			}
		}
		if ($#AlignStart > 10) {	# Arbitrary number to indicate that checking has been succesful
			$checking = 5;
		} elsif ($checking == 0) {	# Let's try the reverse hash and see if we can find matches
			$theHash = \%reverseHash; $genomic = \$QueryRev;
			$checking = 1;
		} else {					# Failure! No matches to sequence
			$checking = 10;
		}
	}
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub WORDCHECK {
	$wordTime = "Find Exons: ";
	$QueryRev = &REVERSE($QuerySeq);
#	&MAKEWORDHASH;
	undef (@Exons); undef (@ExonNames);
	undef (@Introns); undef (@IntronNames);
	@mRNAnames = @Names;
	$ti = time;
	print "<B><FONT COLOR=green>Exon hunting in query with word size of $wordSize. Only the global best matches will be displayed.</FONT></B><BR>\n";
	for $j (0..$#Names) {

		###########################################################################################
		# Search for exact matches of all words in the mRNA (library) with the query (large genomic sequence)
		$checking = 0; $genomic = \$QuerySeq; $theFrame = "Forward";
		while ($checking < 2) {
			undef (@AlignStart); undef (@AlignEnd); undef (@QueryStart); undef (@QueryEnd);
			$RepetitiveWords = 0;
			for $k (0..(length($Sequences[$j]) - $wordSize)) {
				$a = index($$genomic, substr($Sequences[$j], $k, $wordSize));
				if ($a > -1) {
					$b = index($$genomic, substr($Sequences[$j], $k, $wordSize), $a+1 );
					if ($b < 0) {			# Consider the hit only if it is unique (ie, if looking after the match finds no other)
						push @AlignStart,	$k;
						push @AlignEnd,		$k	+ $wordSize-1;
						push @QueryStart,	$a;
						push @QueryEnd,		$a	+ $wordSize-1;
					} else {
						$RepetitiveWords++;
					}
				}
			}
			if ($#AlignStart > length($Sequences[$j])/20) {		# Arbitrary number to indicate that checking has been succesful
				$checking = 5;
			} elsif ($checking == 0) {	# Let's try the reverse hash and see if we can find matches
				$genomic = \$QueryRev;
				$checking = 1; $theFrame = "Reverse";
			} else {					# Failure! No matches to sequence
				$checking = 10;
			}
		}
		
		if ($checking == 10) {
			print "<FONT COLOR=red><B>Only found ",$#AlignStart+1," hits with $Names[$j] using a word size of $wordSize.</B></FONT><BR>\n";
			next;
		} elsif ($RepetitiveWords) {
			$plural = "matches were";
			$plural = "match was" if ($RepetitiveWords == 1);
			print "<FONT COLOR=blue>$RepetitiveWords multi-copy $plural not used.</FONT><BR>\n";
		}
		
		###########################################################################################
		# Collect nearby word matches
		# Rational: Determine the amount of intervening genomic sequence, subtract out the amount of "unclaimed" mRNA sequence.
		# If what is left is less than minIntonSize, assume it should be gathered together into a single exon
		# This allows assembly of regions with less than 100% sequence identity
		$dd = $#AlignStart;	# We will add a dummy record to cap end needed in the while loop below
		print "<FONT COLOR=navy><B>Analyzing ",$dd+1," unique word matches in the <U>$theFrame</U> frame for continuity...</FONT></B><BR>\n";
		$Collected = 0;
		push @AlignStart, 0; push @AlignEnd, 0; push @QueryStart, 999999999; push @QueryEnd, 0;
		for ($k = 1; $k <= $dd; $k++) {									# Gather together matches that are less than an intron away
			$check = $k - 1;
			if ( $QueryStart[$k] - $QueryStart[$check] < -10 ) {		# Allow some negative overlap for small deletions, otherwise disregard as a spurious hit
				print "<FONT COLOR=blue>Ignoring hit to $QueryStart[$k] that preceeds the previous match at $QueryStart[$check].</FONT><BR>\n";
				$AlignStart[$k] = -20;
				next;
			}
			$OkMatch = $check;											# Used for while-loop checking of out-of-order matches
			#      (   intervening genomic sequence    ) - (        Intervening mRNA           )
			while (($QueryStart[$k] - $QueryEnd[$check]) - ($AlignStart[$k] - $AlignEnd[$check]) < $minIntonSize) {	
				if ( $QueryStart[$k] - $QueryStart[$OkMatch] < -10 ) {	# Allow some negative overlap for small deletions, otherwise disregard as a spurious hit
					print "<FONT COLOR=blue>Ignoring hit to $QueryStart[$k] that preceeds the previous match at $QueryStart[$OkMatch].</FONT><BR>\n";
					$AlignStart[$k] = -20;
					$k++;
					next;
				}
				$OkMatch = $k;
				$AlignEnd[$check] = $AlignEnd[$k];
				$QueryEnd[$check] = $QueryEnd[$k];
				$AlignStart[$k] = -20;
				$k++;
			}
			$Collected++;
		}
		
		###########################################################################################
		# Verify reasonable length of final matches, transfer them to a new array
		undef (@Aligned);
		undef (@ExStart); undef (@ExEnd); undef (@QStart); undef (@QEnd);
		print "<FONT COLOR=navy><B>Verifying length and position of $Collected continuous hits...</FONT></B><BR>\n";
		for $k (0..$dd) {
			next if ($AlignStart[$k] < 0);								# This alignment has been merged with a preceeding one.
			if ($AlignEnd[$k] - $AlignStart[$k] + 1 < $minExonSize) {	# The "exon" was probably a single spurious hit and does not represent a real exon
				print "<FONT COLOR=blue>Ignoring likely spurious hit of size ", $AlignEnd[$k] - $AlignStart[$k] + 1, "bp ";
				print "at mRNA position $AlignStart[$k] and query position $QueryStart[$k]</FONT><BR>\n";
				next;
			}
			if ($#Aligned > -1 ) {							# Previous alignments have been analyzed
				($AlignStart, $AlignEnd, $QueryStart, $QueryEnd) = split " ", $Aligned[$#Aligned];
				if ($AlignStart[$k] <= $AlignEnd) {			# Prevent overlap of alignments ( = exons)
					$d = $AlignEnd - $AlignStart[$k] + 1;	# Find the number of bases that are overlapped
					$AlignStart[$k] += $d;					# Lop them off the front of the current alignment
					$QueryStart[$k] += $d;
				}
			}
			push @Aligned, sprintf ("%10.10d %10.10d %10.10d %10.10d",
							$AlignStart[$k], $AlignEnd[$k], $QueryStart[$k], $QueryEnd[$k]);
			push @ExStart,	$AlignStart[$k];				# These are the arrays for the deduced exons
			push @ExEnd,	$AlignEnd[$k];
			push @QStart,	$QueryStart[$k];
			push @QEnd,		$QueryEnd[$k];
		}

		###########################################################################################
		# Check the first and last exons, see if it reasonable to add any additional bases to them
		if ($ExStart[0] < $minExonSize) {						# Remainder to left of first match is too small to be a stand-alone exon
			if ($ExStart[0] > 0) {
				$QStart[0] -= $ExStart[0];						# We'll just tack it on to the first match
				print "<FONT COLOR=blue>Small ($ExStart[0] bp) unmatched 5' fragment appended to first exon.</FONT><BR>\n";
				$ExStart[0] = 0;
			} elsif ($ExStart[0] < 0) {
				print "<FONT COLOR=red><B>WARNING:</B> First exon appears improper: mRNA positions $ExStart[0]-$ExEnd[0], query position $QStart[0]-$QEnd[0].</FONT>\n";
			}
		}
		$siz = length($Sequences[$j]) - 1 - $ExEnd[$#ExEnd];
#		for $d (0..$#ExEnd) {
#			print "$d: mRNA positions $ExStart[$d]-$ExEnd[$d], query position $QStart[$d]-$QEnd[$d].<BR>\n";
#		}
		if ($siz < $minExonSize) {								# Remainder to right of last match is too small to be a stand-alone exon
			if ($siz > 0) {
				$QEnd[$#ExEnd] += $siz;							# Again, we'll just assume it should go with this match
				$ExEnd[$#ExEnd] = length($Sequences[$j])-1;
				print "<FONT COLOR=blue>Small ($siz bp) unmatched 3' fragment appended to last exon.</FONT><BR>\n";
			} elsif ($siz < 0) {
				print "<FONT COLOR=red><B>WARNING:</B> Terminal exon appears improper: mRNA positions $ExStart[$#ExEnd]-$ExEnd[$#ExEnd], query position $QStart[$#ExEnd]-$QEnd[$#ExEnd].</FONT>\n";
			}
		}
		
		###########################################################################################
		# If needed, swap bases between exons to generate better splice sites
		print "<FONT COLOR=navy><B>Optimizing splice junctions for ",$#ExStart + 1," exons...</FONT></B><BR>\n" if ($#ExStart > 0);
		for $k (1..$#ExStart) {					# Check to see if swapping terminal bases on exons results in better splice signals
			$intronL =	$QEnd[$k-1]+1;			# String index of left side of intron
			$intronR =	$QStart[$k]-1;			# String index of right side of intron
			$exonL =	$ExEnd[$k-1];			# Last base on left exon
			$exonR =	$ExStart[$k];			# First base on right exon
			
			$leftLook = 0; $rightLook = 0;
			if ($exonR - $exonL > 1) {			# There is orphaned sequence between the two exons
				$orL = $exonR - $exonL -1;
				$orphan = substr($Sequences[$j],$exonL+1,$orL);
				if ($orL == 1) {				# A single base to deal with
					$ExEnd[$k-1]	+= $orL;	# We'll tack the puppy onto the left exon...
					$exonL			+= $orL;
					$QEnd[$k-1]		+= $orL;
					$intronL		+= $orL;
					$leftLook = -1;				# ... and let the program be allowed to swap it back and forth
				} else {
					if ( substr($orphan, 1) eq substr ($$genomic, $intronL+1, $orL-1)) {	# Is the orphan due to a 5' mismatch on the left exon?
						$ExEnd[$k-1]	+= $orL;	# Extend the left exon by the length of the orphan
						$exonL			+= $orL;
						$QEnd[$k-1]		+= $orL;
						$intronL		+= $orL;
					} elsif ( substr($orphan, 0, $orL-1) eq substr ($$genomic, $intronR-$orL+1, $orL-1)) {
						$ExStart[$k]	-= $orL;	# Extend the right exon by the length of the orphan
						$exonR			-= $orL;
						$QStart[$k]		-= $orL;
						$intronR		-= $orL;
					}
				}
			}
			
			# Can we take bases from the left exon and tack them on the right and still match the underlying genomic sequence?
			while ( substr($Sequences[$j],$exonL+$leftLook,1)	eq substr($$genomic, $intronR+$leftLook, 1) ) { $leftLook--; }
			# As above, but now yanking off the right exon and trying to match them along the left
			while ( substr($Sequences[$j],$exonR+$rightLook,1)	eq substr($$genomic, $intronL+$rightLook,1) ) { $rightLook++; }
			undef (@hits);
			for $j ($leftLook..$rightLook) {
				$score = 2;		# Lower scores are better
				$score-- if ( substr($$genomic, $intronL+$j,2)	=~ m/GT/i);		# Does this arrangement of bases have a nice 5' splice site?
				$score-- if ( substr($$genomic, $intronR+$j-1,2) =~ m/AG/i);	# What about 3'?
				push @hits, sprintf ("%1.1d\t%2.2d\t%d",  $score, abs($j), $j);	# Store the information in a sortable format
			}
			@hits = sort {$a cmp $b} @hits;
			($score,$d,$offset) = split "\t", $hits[0];							# The "lowest" sort is used as the prefered one. It is possible that other hits are equally scored
			if ($ShowIntrons) {
				push @Introns, substr($$genomic, $intronL+$offset, $intronR-$intronL +1);	# Determine the intronic sequence.
				push @IntronNames, ">Intron $k";
			} else {
				$mod = "";
				if ($score > 0) {
					$mod = "<FONT COLOR=brick>Intron $k Sub-ideal: ";
					$mod .= substr($$genomic, $intronL+$offset, 2) . "..." . substr($$genomic, $intronR+$offset-1,2);
					$mod .= "</FONT>";
				}
				$ExonMod[$k-1] = $mod;
			}
			$ExEnd[$k-1] += $offset;
			$ExStart[$k] += $offset;
			$QEnd[$k-1] += $offset;
			$QStart[$k] += $offset;
		}
		undef (@tempSeq); undef (@tempName);
		if ($ExStart[0] > 0) {			# Check for unclaimed sequence at start of mRNA
			push @tempSeq, substr($Sequences[$j], 0, $ExStart[0]);
			push @tempName, ">Orphan 5' Seq";
		}
		for $k (0..$#ExStart) {
			push @tempSeq, substr($Sequences[$j], $ExStart[$k], $ExEnd[$k]-$ExStart[$k] +1);
			push @tempName, ">Exon " . ($k+1);
			if ($k < $#ExStart) {			# Check for unclaimed sequence in middle of mRNA
				if ($ExEnd[$k] + 1 < $ExStart[$k+1]) {
					$orphanCount++;
					push @tempSeq, substr($Sequences[$j], $ExEnd[$k]+1, $ExStart[$k+1] - $ExEnd[$k] -1);
					push @tempName, ">Orphan Internal Seq $orphanCount";
				}
			}
		}
		if ($ExEnd[$#ExStart] < length($Sequences[$j])-1) {			# Check for unclaimed sequence at end of mRNA
			push @tempSeq, substr($Sequences[$j], $ExEnd[$#ExStart]+1);
			push @tempName, ">Orphan 3' Seq";
		}
		push @ExonCount, $#Exons + 1;
		push @Exons, @tempSeq;
		push @ExonNames, @tempName;
	}

	$ti = time - $ti;
	$wordTime .= "$ti sec to assign exon structure.";
	@Sequences = @Exons;
	@Names	= @ExonNames;
	&WRITELIB;
	undef (%forwardHash); undef (%reverseHash);
	undef (@AlignStart); undef (@AlignEnd); undef (@QueryStart); undef (@QueryEnd);
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub COMPARELIBRARY {
	$CompTime = time;
	undef (@Scores); $comparisons = 0;
	# Main loop through library sequences ($i can be considered the "query")
	for $i (0..($#Names-1)) {
		$k = length ($Sequences[$i]) - $libSize;
		next if ($k < 1);
		# Make a hash table for $i
		undef (%LibWord);
		for $n (0..$k) {
			$LibWord{ substr($Sequences[$i], $n, $libSize) } = $n;
		}
		# Loop through all other library sequences ($j's represent the "library" to each "query" $i)
		for $j (($i+1)..$#Names) {
			$k = length ($Sequences[$j]) - $libSize;
			next if ($k < 1);
			$comparisons++;
			# counter will increment with each match (reward matches of increasing size), and will be halved with each non-match (penalize gaps)
			# Score will increase by counter with each match
			$counter = $Score = $RCcounter = $RCScore = 0;
			# Check the hashtable against each word in the forward and reverse complement of $j
			for $n (0..$k) {
				if ($LibWord{ substr($Sequences[$j], $n, $libSize) } ne "") {
					$Score += ++$counter;
				} else {
					$counter /= 2;
				}
				if ($LibWord{ substr($Reverse[$j], $n, $libSize) } ne "") {
					$RCScore += ++$RCcounter;
				} else {
					$RCcounter /= 2;
				}
			}
			$Adjust = 200 / (length($Sequences[$i]) + length($Sequences[$j]));	# Standardize scores to comparison between two 100 bp fragments
			push @Scores, sprintf ("%10.d\t$i\t$j",  int($Score * $Adjust) );
			push @Scores, sprintf ("%10.d\t$i\t$j\tRC",  int($RCScore * $Adjust));
		}
	}
	@Scores = sort {$b cmp $a} @Scores;
	print "<TABLE><TR><TH>Match<TH>Score\n";
	$aCounter = 0;
	for $i (0..$#Scores) {
		($score, $a, $b, $r) = split "\t", $Scores[$i];
		last if ($score < 25 && $i > 30);	# Need to determine what a reasonable cut-off is
		
		$Q = $Sequences[$a]; $L = $Sequences[$b];
		if (length($Sequences[$a]) < length($Sequences[$b])) {
			$c = $a; $a = $b; $b = $c;
		}
		$n = "Q" . $aCounter . "$Here.aln";
		open (FILE, ">$HOME/Data/Alignments/$n") or die "sub COMPARELIBRARY: Failure to write to $n: $!\n";
			print FILE "$Names[$a]\n";
			print FILE "$Sequences[$a]\n";
		close FILE;
		$n = "L" . $aCounter . "$Here.aln";
		open (FILE, ">$HOME/Data/Alignments/$n") or die "sub COMPARELIBRARY: Failure to write to $n: $!\n";
			print FILE "$Names[$b]\n";
			print FILE "$Sequences[$b]\n";
		close FILE;

		$a = substr ($Names[$a],1); $b = substr ($Names[$b],1);		# Chop off ">" from front
		print "<TR><TD>$a <FONT COLOR=blue>vs.</FONT> $b <FONT COLOR=brick>$r</FONT>\n";
		
		print qq(<TD><A HREF="Fasta.cgi?AlignThem=$aCounter&width=500" target=_blank>$score</A>\n);
		$aCounter++;
	}
	print "</TABLE>\n";
	$CompTime = time - $CompTime;
	print "<i>$CompTime sec to perform $comparisons comparisons.</i><BR>\n";
}
# -+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
sub DEFINECONSTANTS {
	$HOME = "/usr/local/etc/httpd/cgi-bin/pagelab/charles";		# The path for Fasta.cgi
	$Here = "."."$ENV{'REMOTE_ADDR'}";							# Returns the IP address of the *USERS* Machine (to allow multiple users)
	$ResLoc =  "/usr/local/etc/httpd/cgi-bin/pagelab/GIF";		# Location for storing gif files
	$ResWeb =  "http://staffa.wi.mit.edu/cgi-bin/pagelab/GIF";											# Where the web should locate the gif files
	
	@RandNucs = ("A","C","T","G",);
	$Allowed = "AGCTUYRWSKMBDHVXN";
	$Reverse = "TCGAARYWSMKVHDBXN";
	%rc = ();
	for $i (0..length($Allowed)-1) {
		$rc{ substr( $Allowed , $i ,1) } = substr( $Reverse , $i ,1 );
	}
	$BlockSize = 80;	# Number of nucs to print in a FastA block
	$hT = (gdTinyFont->height);
	$wT = (gdTinyFont->width);
	$minExonSize = 30;
	$minIntonSize = 30;
	$BigGap = "";
	for $i (1..$minIntonSize) {
		$BigGap .= "-";
	}
	
	$query = new CGI;
	$filename =	$query->param('uploaded_file');
	$libname =	$query->param('uploaded_lib');
	$othername = $query->param('other_file');
	
	@{$Rank{'Zscore'}}	= (0,	30,	50,	100,	200,	300,	350,	400,	450,	500);
	@{$Rank{'SWscore'}}	= (0,	30,	50,	100,	200,	300,	350,	400,	450,	500);
	@{$Rank{'PerID'}}	= (0,	50,	60,	70,		80,		85,		90,		95,		98,		100);
#	@{$Rank{'expect'}}	= (2,	0,	-2,	-4,		-8,		-12,	-14,	-18,	-25,	-30);
	
	@dataTypes = (	"Name",		"initn",	"init1",	"opt",			"Zscore",	"expect",
					"SWscore",	"PerID",	"Overlap",	"MatchStart",	"MatchEnd",	"LibStart",
					"LibEnd",	"Files",	);
	
	# The following routine retrieves values from the CGI paramater lists
	# It also sets defaults if no value is specified.				
	@CGIParams = (	"MinPalSize=4",		"MaxPalSize=60",	"MinSimple=9",		"SimpStep=1",	"SimpSize=6",
					"MinLoop=4",		"MaxLoop=10",		"MinHairpin=10",	"BinSize=50",
					"DoCPG=0",			"DoPalindrome=0",	"ChunkSize=20000",	"ChunkOverlap=40",
					"minBinSize=250",	"binFlank=1",		"FourFrequency=",	"width=750",
					"CheckRepeat=",		"SortBy=SWscore",	"FindExons=",		"ShowIntrons=",
					"wordSize=15",		"Cutoff=1",			"LibCompare=",		"libSize=8",
					"ColorType=",					);
	foreach $i (@CGIParams) {	
		($var, $val) = split "=", $i;
		$$var = $query->param($var) if ($query->param($var) ne "");
		$$var = $val if ($$var eq "");
	}
	
	#	MinPalSize	=	Maximum reach to consider
	#	MaxPalSize	= 	Minimum reach to consider 0 = 2bp palindromes
	#	MinSimple	=	Minimum size to be considered a simple repeat
	#	SimpStep	=	Steps to take while looking for simple repeats (suggest 1 or 2)
	#	SimpSize	=	Maximum unit size of simple repeat to consider
	#	MinLoop		=	Minimum number of nucleotides to consider in the loop of a hairpin
	#	MaxLoop		=	Maximum number of nucleotides to consider in the loop of a hairpin
	#	MinHairpin	=	Minimum number of bases needed to be considered a hairpin
	#	BinSize		=	Size of window to analyze palindrome frequency
	#	ChunkSize	=	Maximum fragment size to pass to FastA
	$ChunkSize	= 20000 if ($ChunkSize > 20000);	# This is the maximum currently allowed by FastA
	$MinReach = $MinPalSize/2 - 1;
	$MaxReach = $MaxPalSize/2 - 1;
	
	# Minimum number of repeats to be considered valid
	@RepeatNames =	("Null", "Simple", "Dinucleotide", "Trinucleotide", "Tetranucleotide", "Pentanucleotide", "Hexanucleotide",	"Heptanucleotide");
	@SimpLimit =	(10000,		8,			5,				4,					3,				3,					3,					3);
	@F = (); 
	print "<!"; system "rm $ResLoc/*graph*$Here.tmp"; print ">";	# Remove old gif files
	$citation = "FASTA citation: W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448";

}
