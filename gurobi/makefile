# https://support.gurobi.com/hc/en-us/articles/360039093112
# Gurobi MakeFile for C++

PLATFORM = linux64
INC      = /home/morete/ic/gurobi903/include
CPP      = g++ -std=c++17
CARGS    = -m64 -g
CPPLIB   = -L/home/morete/ic/gurobi903/lib -lgurobi_c++ -lgurobi90
GRBAPP   = DOTNETCore2

STD = -std=c++17
OP1 = -Wall -Wextra -Wpedantic -Wshadow
OP2 = -Wno-unused-result -Wno-unused-function
OP3 = -Wfloat-equal -Wconversion -Wshift-overflow=2 -Wlogical-op
SAN = -fsanitize=address -fsanitize=undefined -fno-sanitize-recover

%: %.cpp
	$(CPP) $(CARGS) $(SAN) -o $@ $< -I$(INC) $(CPPLIB) -lm

clean:
	rm -rf *.o *_c *_c++ *.class *.log *.rlp *.lp *.bas *.ilp *.mps *.prm; \
	if [ -d $(GRBAPP) ]; then \
		cd $(GRBAPP); \
		find . ! \( -name dotnetcore2.csproj -o -name . \) -exec rm -rf {} +; \
	fi
