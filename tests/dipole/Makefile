all:
	@make tests/dipole_executables -C ../..

r: results
results:
	@make tests/dipole_clean_results -C ../..

c: clean
clean: results
	@make tests/dipole_clean -C ../..
