#!/usr/bin/perl -w

require LWP::UserAgent;

my $ua = LWP::UserAgent->new;
$ua->timeout(60);

sub toplist {
  my $response = $ua->get("http://csclub.uwaterloo.ca/contest/rankings.php");
  my $c = $response->content;
  my @users;
  while($c =~ /user_id=(\d+)(.*)$/ms) {
    $c = $2;
    push @users, $1;
  }
  return @users;
}

sub gamelist {
  my ($userid) = @_;
  my $response = $ua->get("http://csclub.uwaterloo.ca/contest/profile_games.php?user_id=$userid");
  my $c = $response->content;
  my @games;
  while($c =~ /game_id=(\d+)(.*)$/ms) {
    $c = $2;
    push @games, $1;
  }
  return @games;
}

my @users = toplist;
print("# loaded " . scalar(@users). " users\n");

map { map { print("./getgame.pl $_; echo $_ >>allgames\n") } gamelist($_) } @users;

