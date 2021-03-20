.PHONY: clean All

All:
	@echo "----------Building project:[ mirseek - Release ]----------"
	@cd "mirseek" && "$(MAKE)" -f  "mirseek.mk"
clean:
	@echo "----------Cleaning project:[ mirseek - Release ]----------"
	@cd "mirseek" && "$(MAKE)" -f  "mirseek.mk" clean
