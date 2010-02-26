#!/usr/bin/perl -w
use List::Util qw(sum);
# model: v1-v0 = A [n0 e0 a0 n1 e1 a1]

my ($denom,$n) = (0,0);
my @num = (0,0,0);
while(<>) {
  my ($v0, $v1, $st, $n0,$e0,$a0, $n1,$e1,$a1) = split / /;
  pop @params; # remove redundant distance parameter
  if($n0+$e0+$a0 == $n1+$e1+$a1) { next; }
  my $y = $v1-$v0;
  @params = ($n1-$n0, $e1-$e0, $a1-$a0);
  $denom += sum(map { $_*$_ } @params);
  for($i=0;$i<@num;$i++) { $num[$i] += $params[$i] * $y; }
  $n++;
}

print("$n datapoints; denom = $denom\n");
for($i=0;$i<@num;$i++) { print 1000*$num[$i]/$denom . " "; }
print("\n");

