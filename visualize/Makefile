.PHONY: all
all: out/index.js

out/index.js: src/index.ts build.ts
	deno run --allow-env --allow-read --allow-write build.ts
