#!/usr/bin/perl -w

require LWP::UserAgent;

my $ua = LWP::UserAgent->new;
$ua->timeout(10);

if(!defined $ARGV[0]) {
  die "usage: $0 <game_id>";
}

my $flip = $ARGV[1];
my $timeorder = $ARGV[2];

my $response = $ua->get("http://csclub.uwaterloo.ca/contest/generate_game_summary.php?game_id=$ARGV[0]");

my ($ok,$mapsize,$map,$p1,$p2,$moves,$ok2) = split(/\|/, $response->content);
my @moves = split(//, $moves);

my ($w,$h) = split(/\s+/, $mapsize);
my %map;
my $idx=0;
my ($p1x, $p1y, $p2x, $p2y);
$map = join('', split(/\n/, $map));
for(my $j=0;$j<$h;$j++) {
  for(my $i=0;$i<$w;$i++,$idx++) {
    my $c = substr($map, $idx, 1);
    if($c eq "1") { $p1x = $i; $p1y = $j; }
    if($c eq "2") { $p2x = $i; $p2y = $j; }
    $map{$i,$j} = substr($map, $idx, 1);
  }
}

sub _next {
  my ($x, $y, $move) = @_;
  if($move eq 'N') { return ($x, $y-1); }
  if($move eq 'S') { return ($x, $y+1); }
  if($move eq 'E') { return ($x+1, $y); }
  if($move eq 'W') { return ($x-1, $y); }
}

sub move {
  my ($m1, $m2) = @_;
  $map{$p1x,$p1y} = '#';
  $map{$p2x,$p2y} = '#';
  ($p1x,$p1y) = _next($p1x,$p1y,$m1);
  ($p2x,$p2y) = _next($p2x,$p2y,$m2);
  $map{$p1x,$p1y} = '1';
  $map{$p2x,$p2y} = '2';
}

sub display_map {
  my @out = ("$mapsize\n");
  for(my $j=0;$j<$h;$j++) {
    for(my $i=0;$i<$w;$i++,$idx++) {
      my $c = $map{$i,$j};
      if($flip && $c eq '1') { $c = '2' }
      elsif($flip && $c eq '2') { $c = '1' }
      push @out, $c;
    }
    push @out, "\n";
  }
  return join('', @out);
}

#print(" --- $p1 ($p1x,$p1y) vs $p2 ($p2x,$p2y)\n");
my @out = (display_map);

while(@moves) {
  my $m1 = shift @moves;
  my $m2 = shift @moves;
  move($m1,$m2);
  last if($p1x == $p2x && $p1y == $p2y);
  push @out, display_map();
}
pop @out; # last move is often confusing my examiner

if($timeorder) { @out = reverse @out; }

print join('', @out);

