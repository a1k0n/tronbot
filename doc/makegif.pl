use strict;
use GD;
use Data::Dumper;

my ($w,$h) = (16,15);
my $scale = 8;

# setup the image
my $image = GD::Image->new($w*$scale,$h*$scale);
my $white = $image->colorAllocate(255,255,255);
my $black = $image->colorAllocate(0,0,0);

# thanks, adobe kuler
my $player1 = $image->colorAllocate(0,145,229);
my $area1 = $image->colorAllocate(169,223,255);
my $player2 = $image->colorAllocate(229,0,7);
my $area2 = $image->colorAllocate(255,188,177);

my %colormap = (
  "#" => $black,
  "A" => $player1,
  "1" => $area1,
  "o" => $area1,
  "B" => $player2,
  "2" => $area2,
  "x" => $area2,
  "." => $white
);

#$image->transparent($white);

# setup some font goodies
my $fontcolor = $image->colorAllocate(0,0,0);
my $font = GD::Font->Small();

# setup some settings into variables
my $loop = 1;
my $speed = 100; # 1/100 of a sec
my $x_font = 40; # from right (x or y ??)
my $y_font = 40; # from top (x or $y ??)

my @frames = ();
while(<>) {
  if(/evaluating board/) { push @frames, []; }
  elsif(s/^~~~ //) {
    chomp;
    my @line = split(//);
    push @{$frames[@frames-1]}, \@line;
  }
}

sub frame_handler {
  my ($im, $frame) = @_;
  for(my $j=0;$j<@$frame;$j++) {
    my $line = $frame->[$j];
    for(my $i=0;$i<@$line;$i++) {
      my $ch = $line->[$i];
      my $c = $colormap{$ch};
      $im->filledRectangle($scale*$i, $scale*$j, $scale*($i+1)-1, $scale*($j+1)-1, $c);
      if($ch eq 'x' || $ch eq 'o') {
        print "cut vertex @ $i $j\n";
        $im->arc($scale*($i+0.5), $scale*($j+0.5), $scale*0.75, $scale*0.75, 0, 360, $black);
      }
    }
  }
}

my $n=0;
foreach my $framedata (@frames) {
   # make a frame of right size
   $n++;
   open GIF, sprintf(">board%03d.png", $n);
   binmode GIF;
   frame_handler($image, $framedata);              # add the data for this frame
   print GIF $image->png;
}

