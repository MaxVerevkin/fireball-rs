import { transpile } from "https://deno.land/x/emit@0.40.0/mod.ts";

const url = new URL("./src/index.ts", import.meta.url);
const result = await transpile(url);

const code = result.get(url.href);
// console.log(code)
await Deno.writeTextFile("./out/index.js", code);
