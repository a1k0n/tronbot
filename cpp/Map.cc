// Map.cc

#include "Map.h"
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

Map::Map() {
  ReadFromFile(stdin);
}

int Map::Width() const {
  return map_width;
}

int Map::Height()  const {
  return map_height;
}

bool Map::IsWall(int x, int y) const {
  if (x < 0 || y < 0 || x >= map_width || y >= map_height) {
    return true;
  } else {
    return is_wall[x][y];
  }
}

int Map::MyX() const {
  return player_one_x;
}

int Map::MyY() const {
  return player_one_y;
}

int Map::OpponentX() const {
  return player_two_x;
}

int Map::OpponentY() const {
  return player_two_y;
}

void Map::MakeMove(const std::string& move) {
  if (move.length() == 0) {
    fprintf(stderr, "ERROR: zero-length string passed to MakeMove(string)\n");
    MakeMove(0);
  } else {
    int c = (int)move[0];
    switch (c) {
    case 'n':
    case 'N':
      MakeMove(1);
      break;
    case 'e':
    case 'E':
      MakeMove(2);
      break;
    case 's':
    case 'S':
      MakeMove(3);
      break;
    case 'w':
    case 'W':
      MakeMove(4);
      break;
    default:
      fprintf(stderr, "Invalid string passed to MakeMove(string): %s\n"
	      "Move string must start with N, E, S, or W!", move.c_str());
      MakeMove(0);
      break;
    }
  }
}

void Map::MakeMove(int move) {
  fprintf(stdout, "%d\n", move);
  fflush(stdout);
}

void Map::ReadFromFile(FILE *file_handle) {
  int x, y, c;
  int num_items = fscanf(file_handle, "%d %d\n", &map_width, &map_height);
  if (feof(file_handle) || num_items < 2) {
    exit(0); // End of stream means end of game. Just exit.
  }
  is_wall =
    std::vector<std::vector<bool> >(map_width,
				    std::vector<bool>(map_height, false));
  x = 0;
  y = 0;
  while (y < map_height && (c = fgetc(file_handle)) != EOF) {
    switch (c) {
    case '\r':
      break;
    case '\n':
      if (x != map_width) {
	fprintf(stderr, "x != width in Board_ReadFromStream\n");
	return;
      }
      ++y;
      x = 0;
      break;
    case '#':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      is_wall[x][y] = true;
      ++x;
      break;
    case ' ':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      is_wall[x][y] = false;
      ++x;
      break;
    case '1':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      is_wall[x][y] = false;
      player_one_x = x;
      player_one_y = y;
      ++x;
      break;
    case '2':
      if (x >= map_width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      is_wall[x][y] = false;
      player_two_x = x;
      player_two_y = y;
      ++x;
      break;
    default:
      fprintf(stderr, "unexpected character %d in Board_ReadFromStream", c);
      return;
    }
  }
}
