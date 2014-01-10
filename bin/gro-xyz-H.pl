#!/usr/bin/perl

$i = 0;
while(<>){
	if(/^ *[0-9]+([A-Z][A-Z][A-Z]*) +(.+) +(\-*[0-9.]+) +(\-*[0-9.]+) +(\-*[0-9.]+)\n/){
		$a1 = $1;
		$a2 = $2;
		$a3 = $3;
		$a4 = $4;
		$a5 = $5;
		if ($a1 =~ /SOL/) {$atom_temp[$i] = 10;}
		else {$atom_temp[$i] = 1;}
		$atom_x[$i] = $a3;
		$atom_y[$i] = $a4;
		$atom_z[$i] = $a5;
		$atom_wt[$i] = 1;
		if($a2=~ m/^C/){$atom_num[$i] = 12;}
		if($a2=~ m/^O/){$atom_num[$i] = 16;}
		if($a2=~ m/^N/){$atom_num[$i] = 14;}
		if($a2=~ m/^P/){$atom_num[$i] = 31;}
		if($a2=~ m/^S/){$atom_num[$i] = 32;}
		if($a2=~ m/^H/ | $a2=~ m/^'H/){$atom_num[$i] = 1;}
		$i++;
	}
	if(/^  [0-9.]+ +[0-9.]+ +[0-9.]+/){
		split;
		$x = $_[0] * 10;
		$y = $_[1] * 10;
		$z = $_[2] * 10;
	} 
}
print "Title\n";
print $x." ".$y." ".$z."\n";
for($c = 0; $c < $i; $c++){
		$x = $atom_x[$c] *10;
		$y = $atom_y[$c] *10;
		$z = $atom_z[$c] *10;
		printf("%d %g %g %g %g %g\n", $atom_num[$c],$x,$y,$z,$atom_wt[$c],$atom_temp[$c],);
}
print "-1\n";
