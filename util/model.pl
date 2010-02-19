#!/usr/bin/perl -w
use List::Util qw(sum);
# model: (v1-v0)/st = A [n0 e0 a0 n1 e1 a1 log(allpaths)]
# e.g.
# v1 v0 st n0 e0  a0 n1 e1  a1 d  np
# 54 50 14 56 185 6  46 151 10 12 152
# 54 50 15 85 291 7  48 155 11 10 55
# 54 50 16 79 270 6  56 186 11 10 44
# 54 50 17 84 289 8  53 176  9 12 180
# 54 50 18 83 286 7  56 188  9 14 1318
# 54 50 19 72 248 6  63 215 10 16 2380
# 54 50 20 75 260 7  67 230 10 18 3807
# 54 50 21 72 249 9  72 249  9 20 11460
# 54 50 22 73 252 8  73 252  8 22 11460
# 54 50 23 74 256 10 74 256 10 24 129828

my ($denom,$n) = (0,0);
my @num = (0,0);
while(<>) {
  my ($v0, $v1, $st, $n0,$e0,$a0, $n1,$e1,$a1, $d, $np) = split / /;
  pop @params; # remove redundant distance parameter
  if($n0+$e0+$a0 == $n1+$e1+$a1) { next; }
  my $lnp = log($np+1);
  my $y = ($v1-$v0)/(1+$st);
  @params = (($n1-$n0)/$lnp, ($e1-$e0)/$lnp);
  $denom += sum(map { $_*$_ } @params);
  for($i=0;$i<@num;$i++) { $num[$i] += $params[$i] * $y; }
  $n++;
}

print("$n datapoints; denom = $denom\n");
for($i=0;$i<@num;$i++) { print 1000*$num[$i]/$denom . " "; }
print("\n");

