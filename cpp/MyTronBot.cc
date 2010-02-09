// MyTronBot.cc
//
// This is the code file that you will modify in order to create your entry.

#include "Map.h"
#include <string>
#include <vector>

std::string MakeMove(const Map& map) {
  int x = map.MyX();
  int y = map.MyY();
  if (!map.IsWall(x, y-1)) {
    return "NORTH";
  }
  if (!map.IsWall(x+1, y)) {
    return "EAST";
  }
  if (!map.IsWall(x, y+1)) {
    return "SOUTH";
  }
  if (!map.IsWall(x-1, y)) {
    return "WEST";
  }
  return "NORTH";
}

// Ignore this function. It is just handling boring stuff for you, like
// communicating with the Tron tournament engine.
int main() {
  while (true) {
    Map map;
    Map::MakeMove(MakeMove(map));
  }
  return 0;
}
