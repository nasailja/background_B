INCLUDE_FILE = makefiles/macosx_macports

CXXFLAGS = -std=c++0x -W -Wall -Wextra -pedantic -O3

EXECUTABLES = \
  tests/integration/odeint.bexe \
  tests/integration/cubature.exe

TESTS = \
  tests/integration/odeint.btst \
  tests/integration/cubature.tst

include $(INCLUDE_FILE)

%.exe: %.cpp Makefile
	@echo "CXX "$< && $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(CUBATURE_CPPFLAGS) $(CUBATURE_LDFLAGS) $(CUBATURE_LIBS) $< -o $@

# these require boost
%.bexe: %.cpp Makefile
	@echo "CXX "$< && $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(BOOST_CPPFLAGS) $< -o $@


%.tst: %.exe
	@echo "RUN "$< && ./$< && echo PASS && touch $@

%.btst: %.bexe
	@echo "RUN "$< && ./$< && echo PASS && touch $@


all: test

t: test
test: $(EXECUTABLES) $(TESTS)
	@echo && echo "All tests passed."

c: clean
clean:
	@echo "CLEAN" && rm -f $(EXECUTABLES) $(TESTS)
