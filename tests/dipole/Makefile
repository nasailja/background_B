all:
	@make tests/dipole_tests -C ../..

r: results
results:
	@make tests/dipole_clean_results -C ../..

c: clean
clean: results
	@make tests/dipole_clean -C ../..
