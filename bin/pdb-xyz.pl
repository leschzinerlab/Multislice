#!/usr/bin/perl

$i = 0;
while(<>){
	if(/^CRYST1/){
		split;
		$x = $_[1];
		$y = $_[2];
		$z = $_[3];
	}
	if(/^ATOM/ ){
		split;
		$a1 = $_[3];
		$a2 = $_[2];
		$a3 = $_[5];
		$a4 = $_[6];
		$a5 = $_[7];
		if ($a1 =~ /SOL/) {$atom_temp[$i] = 800;}
		else {$atom_temp[$i] = 80;}
		$atom_x[$i] = $a3;
		$atom_y[$i] = $a4;
		$atom_z[$i] = $a5;
		$atom_wt[$i] = 1;
		if($a2=~ m/^C/){$atom_num[$i] = 12;}
		if($a2=~ m/^O/){$atom_num[$i] = 16;}
		if($a2=~ m/^N/){$atom_num[$i] = 14;}
		if($a2=~ m/^P/){$atom_num[$i] = 31;}
		if($a2=~ m/^S/){$atom_num[$i] = 32;}
		if($a2=~ m/^H/){$atom_num[$i] = 1;}
		$i++;
	}
}
print "Title\n";
print $x." ".$y." ".$z."\n";
for($c = 0; $c < $i; $c++){
		$x = $atom_x[$c];
		$y = $atom_y[$c];
		$z = $atom_z[$c];
		printf("%d %g %g %g %g %g\n", $atom_num[$c],$x,$y,$z,$atom_wt[$c],$atom_temp[$c],);
}
print "-1\n";
