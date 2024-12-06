#!/usr/local/bin/perl
BEGIN {
    print "Content-type: text/html\n\n";
    # Ensure that errors go to the web browser. From Steve's Primer3 program
    open(STDERR, ">&STDOUT");
    $| = 1;
    print '';
}

# Staffa path for this file is:
# /usr/local/etc/httpd/cgi-bin/pagelab/charles/
# http://staffa.wi.mit.edu/cgi-bin/pagelab/charles/blank.cgi

strict;
use CGI;
use GD;
print <<EOF;
<html><TITLE>Marker Oligo Number Finder</TITLE>
<BODY BGCOLOR=white  text="#000000" link="#118911">
EOF

require "qd.pl";
require "ParseQuery.pl";
&PARSE;
&LOADDATA;

print $query->start_multipart_form(
							-action=>"blank.cgi",
							-method=>POST,
							-name=>"MainForm");
print "Specify a subset of markers to list here:<BR>";
print $query->textarea(		-name=>'myMarkers',
							-default=>$myMarkers,
							-title=>"You can paste a list of marker names (separated by returns) into this field to analyze just a subset of the vector data.",
							-rows=>12,
							-columns=>25);
print $query->submit(		-value=>'Execute',
							-title=>"Run the program with the settings and files indicated above.",
							),"<BR>\n";

print $query->endform;

@markers = ();
foreach $n (keys %left) {
	next if ($n eq "");
	$_ = $n;
	if (/sY/) {
		s/sY//;
		push @markers, sprintf("%6.6d",$_);
	} else {
		push @markers, $_;
	}
}
@markers = sort {$a cmp $b} @markers;
for $i (0..$#markers) {
	if ($markers[$i] > 0) {
		
		$markers[$i] = sprintf("sY%d" , $markers[$i]+0);
	}
}
if ($myMarkers ne "") {
	@tempOrder = split /[\r\n]/, $myMarkers;
	$targ = \@tempOrder;
} else {
	$targ = \@markers;
}
print "<PRE>\n\nlocus\toligo l\toligo r\tsize\toligo sequence left\toliogo sequence right\n";
foreach $sts (@{$targ}) {
	next if ($sts eq "" || $sts eq "000000");
	print"$sts\t$left{$sts}\t$right{$sts}\t$size{$sts}\t$leftSeq{$sts}\t$rightSeq{$sts}\n";
}
print "</PRE>";
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
sub LOADDATA {
	$HOME = "/usr/local/etc/httpd/cgi-bin/pagelab/charles";
	$HELEN = "/usr/local/etc/httpd/cgi-bin/pagelab";
	$Here = "."."$ENV{'REMOTE_ADDR'}";				# Returns the IP address of the *USER'S* Machine (to allow multiple users)
	my $dataLoc = "~Tomoko/";
	my $priLoc	= $dataLoc . "map/";
	my $secLoc	= $dataLoc . "wu_data/database/";
	$query = new CGI;
	
	open (FILE, "$HELEN/sts_data") or die "Failure to open $HELEN/sts_data: $!\n";
	while (<FILE>) {
		chomp;
		my ($action, $value) = split " ", $_;
		$sts = $value if ($action =~ /STS/i);
		$left{$sts} = $value if ($action =~ /Oligo_left/i);
		$right{$sts} = $value if ($action =~ /Oligo_right/i);
		$leftSeq{$sts} = $value if ($action =~ /Sequence_left/i);
		$rightSeq{$sts} = $value if ($action =~ /Sequence_r/i);
		$size{$sts} = $value if ($action =~ /size/i);
	}
	close FILE;
	
	@CGIParams = ("myMarkers=",);
	foreach $i (@CGIParams) {	
		($var, $val) = split "=", $i;
		$$var = $query->param($var) if ($query->param($var) ne "");
		$$var = $val if ($$var eq "");
		$GV{$var} = $$var;
	}
	
}
# + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - + - 
