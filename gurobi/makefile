# https://support.gurobi.com/hc/en-us/articles/360039093112
# Gurobi MakeFile for C++
# Based on Gurobi exemple

PLATFORM = linux64
INC      = /home/morete/ic/gurobi903/include
CPP      = g++ -std=c++17
CARGS    = -m64 -g
CPPLIB   = -L/home/morete/ic/gurobi903/lib -lgurobi_c++ -lgurobi90
GRBAPP   = DOTNETCore2

%: %.cpp
	$(CPP) $(CARGS) -o $@ $< -I$(INC) $(CPPLIB) -lm

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name dotnetcore2.csproj -o -name . \) -exec rm -rf {} +; \
	fi
