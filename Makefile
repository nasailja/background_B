# to compile and test pick one makefile from makefiles directory
# or write a new one for your environment
ENVIRONMENT_MAKEFILE = makefiles/macosx_macports

#
# The lines below are not intended to be modified by users
#
CXXFLAGS = -std=c++0x -W -Wall -Wextra -pedantic -O3
CPPFLAGS = -I source
include $(ENVIRONMENT_MAKEFILE)

EXECUTABLES =
TESTS =
RESULTS =
CLEAN =

all: test

include \
  tests/integration/project_makefile \
  tests/dipole/project_makefile

%.tst: %.exe
	@echo "RUN "$< && ./$< && echo PASS && touch $@

t: test
test: $(EXECUTABLES) $(TESTS)
	@echo && echo "All tests passed."

r: results
results: $(RESULTS)
	@echo results

c: clean
clean: results $(CLEAN)
	@echo clean
