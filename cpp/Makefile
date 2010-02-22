CXXFLAGS = -g -O3 -funroll-loops -Wall

OBJS = MyTronBot.o
REFBOT = ./MyTronBot8

MyTronBot: $(OBJS)
	$(CXX) $(CXXFLAGS) -o MyTronBot $(OBJS)

clean::
	rm -f $(OBJS) MyTronBot a1k0n.zip

test:: MyTronBot
	java -jar engine/Tron.jar maps/apoc_small.txt ./MyTronBot "java -jar example_bots/Chaser.jar" 0 2>&1

submission: MyTronBot.cc
	zip a1k0n.zip MyTronBot.cc artictbl.h

testall:: MyTronBot
	rm -f testresults.out
	for i in maps/*.txt; do java -jar engine/Tron.jar $$i ./MyTronBot $(REFBOT) 0 2>&1 | tee -a testresults.out; done

testall2:: MyTronBot
	rm -f testresults2.out
	for i in maps/*.txt; do java -jar engine/Tron.jar $$i $(REFBOT) ./MyTronBot 0 2>&1 | tee -a testresults2.out; done

