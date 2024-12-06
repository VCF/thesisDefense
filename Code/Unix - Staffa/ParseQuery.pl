# Query Parse subroutine
# Parses the passed query string and assigns values to
# parameters designated by the string.
# Formats recognized:
#	x,y					-> $x = x, $y = y
#	name = a'b'c'd'...	-> @name = (a,b,c,d,...)
#	name@x = a			-> $name[x] = a
#	name = x			-> $name = x
$| = 1;

sub PARSE {
	$qs = $ENV{'QUERY_STRING'};
	
	$_ = $qs;
	s/\x2B/ /g;		# Space
	s/%20/\x20/g;	# Space
	s/%22/\x22/g;	# "
	s/%27/\x27/g;	# '
	s/%28/\x28/g;	# (
	s/%29/\x29/g;	# )
	s/%2B/\x2B/g;	# +
	s/%2D/\x2D/g;	# -
	s/%2F/\x2F/g;	# /
	s/%3A/\x3A/g;	# :
	s/%3B/\x3B/g;	# ;
	s/%3C/\x3C/g;	# <
	s/%3E/\x3E/g;	# >
	s/%3F/\x3F/g;	# ?
	s/%7C/\x7C/g;	# |
	s/%40/\x40/g;	# @
	s/%5B/\x5B/g;	# [
	s/%5C/\x5C/g;	# ]
	s/%0D%0A/\n/g;	# Return
	s/%09/\t/g;		# tab
	
	$qs = $_;
	
	if ("null" eq ".submit") {	# Data was submitted with CGI standard form
		@inp = split "\n",$qs; chomp (@inp);
		for $i (0..$#inp) {
			$_ = $inp[$i]; s/%2C/\x2C/g; s/%26/\x26/g; $inp[$i] = $_;		#Swaps , and &
			($b, $c) = split "=",$inp[$i];
			next unless ($b);
			$$b = $c;
		}
	} else {
		@inp = split "&",$qs; chomp (@inp);
		for $i (0..$#inp) {
			if ($inp[$i] =~ m/,/) {
				($x, $y) = split ",",$inp[$i];
			} else {
				$_ = $inp[$i]; s/%2C/\x2C/g; s/%26/\x26/g; $inp[$i] = $_;		#Swaps , and &
				($b, $c) = split "=",$inp[$i];
				if ($inp[$i] =~ m/'/) {
					@$b = (split "'", $c);
				} elsif ($b =~ "BOX") {		#BOX designates grouped checkbox arrays in forms
					push @$b, $c;
				} elsif ($b =~ "HASH") {	#HASH designates member of hash
					($d,$e) = split "HASH", $b;
					$$d{$e} = $c;
				} elsif ($b =~ m/\@/ ) {
					($d,$e) = split "\@",$b;
					$$d[$e] = $c;
				} else {
					$$b = $c;
				}
			}
		}
	}
}
