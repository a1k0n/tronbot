CXXFLAGS = -g -O2 -Wall

OBJS = MyTronBot.o

MyTronBot: $(OBJS)
	$(CXX) $(CXXFLAGS) -o MyTronBot $(OBJS)

clean::
	rm -f $(OBJS) MyTronBot a1k0n.zip

test:: MyTronBot
	java -jar engine/Tron.jar maps/apocalyptic.txt ./MyTronBot "java -jar example_bots/WallHugger.jar" 0 2>&1

submission: MyTronBot.cc
	cd ..; zip cpp/a1k0n.zip cpp/MyTronBot.cc cpp/Makefile

testall:: MyTronBot
	for i in maps/*.txt; do java -jar engine/Tron.jar $$i ./MyTronBot "java -jar example_bots/WallHugger.jar" 0 2>&1; done

