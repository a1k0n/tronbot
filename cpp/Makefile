CXXFLAGS = -g -O3 -Wall

OBJS = MyTronBot.o

MyTronBot: $(OBJS)
	$(CXX) $(CXXFLAGS) -o MyTronBot $(OBJS)

clean::
	rm -f $(OBJS) MyTronBot

test:: MyTronBot
	java -jar engine/Tron.jar maps/apocalyptic.txt ./MyTronBot "java -jar example_bots/WallHugger.jar" 0 2>&1

submission: MyTronBot.cc
	cd ..; zip cpp/a1k0n.zip cpp/MyTronBot.cc cpp/Makefile

test2:: MyTronBot
	java -jar engine/Tron.jar maps/u.txt ./MyTronBot "java -jar example_bots/WallHugger.jar" 0 2>&1
