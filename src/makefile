# Gurobi MakeFile for C++
# Based on Gurobi exemple

PLATFORM = linux64
INC      = /home/morete/ic/gurobi/include
CPP      = g++ -std=c++17
CARGS    = -m64 -g
CPPLIB   = -L/home/morete/ic/gurobi/lib -lgurobi_c++ -lgurobi90
GRBAPP   = DOTNETCore2
OPT      = -Ofast
%: %.cpp
	$(CPP) $(CARGS) $(OPT) -o $@ $< -I$(INC) $(CPPLIB) -lm

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name dotnetcore2.csproj -o -name . \) -exec rm -rf {} +; \
	fi
