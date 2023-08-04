# Extended DIMACS format

In the extended DIMACS format, a literal is represented by the string `[+-]idx:value`:

- `+idx:value` or `idx:value` represents the literal _x >= value_, with x being the variable represented by the index `idx`.

- `-idx:value` represents _x <= value_.

Clauses still end with the value "0".

In an extended DIMACS formula, a literal of the form `[+-]idx` without the value specified is interpreted as either _x >= 1_ or _x <= 0_, depending on the sign.
