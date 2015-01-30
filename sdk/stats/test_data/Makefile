CPPFLAGS = -Wno-deprecated -Wall -pedantic
LDFLAGS = -lm -L $(R_LD_LIBRARY_PATH) -lRmath

BASENAME = make_test_data

test_data: $(BASENAME)
	./$(BASENAME)
	R --no-save < make_test_data.R

$(BASENAME): $(BASENAME).o
	$(CXX) -o $(BASENAME) $(BASENAME).o $(LDFLAGS) 

$(BASENAME).o: $(BASENAME).cpp
	$(CXX) $(CPPFLAGS) -I $(R_INCLUDE) -DMATHLIB_STANDALONE -c $(BASENAME).cpp
	
clean:
	rm -rf *.o $(BASENAME)

