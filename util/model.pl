#!/usr/bin/perl -w
use List::Util qw(sum);
# model: (v1-v0) = A [n0 e0 a0 n1 e1 a1 dist]

my ($denom,$n) = (0,0);
my @num = (0,0,0);
while(<>) {
  my ($v0, $v1, $turn, @params) = split / /;
  pop @params; # remove redundant distance parameter
  my $dist = pop @params;
  my ($n0,$e0,$a0, $n1,$e1,$a1) = @params;
  if($n0+$e0+$a0 == $n1+$e1+$a1) { next; }
  @params = ($n1-$n0, $e1-$e0, $a1-$a0);
  $denom += sum(map { $_*$_ } @params);
  for($i=0;$i<@num;$i++) { $num[$i] += $params[$i] * ($v1-$v0); }
  $n++;
}

print("$n datapoints; denom = $denom\n");
for($i=0;$i<@num;$i++) { print 1000*$num[$i]/$denom . " "; }
print("\n");

