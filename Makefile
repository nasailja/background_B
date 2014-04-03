include makefiles/macosx_macports

CXXFLAGS = -std=c++0x -W -Wall -Wextra -pedantic -O3

EXECUTABLES = \
  tests/integration/odeint.bexe

TESTS = \
  tests/integration/odeint.btst


# these require boost
%.bexe: %.cpp Makefile
	@echo "CXX "$< && $(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(BOOST_CPPFLAGS) $< -o $@


%.btst: %.bexe
	@echo "RUN "$< && ./$< && echo PASS && touch $@


all: test

t: test
test: $(EXECUTABLES) $(TESTS)
	@echo && echo "All tests passed."

c: clean
clean:
	@echo "CLEAN" && rm -f $(EXECUTABLES) $(TESTS)
