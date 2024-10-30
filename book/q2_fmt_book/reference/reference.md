# q2-FMT Commands

We are currently working on a way to auto generate help text and examples associated with q2-fmt.
For now the best way to get information about which commands are available in q2-fmt is to run `qiime2 fmt --help`.
For more information on each command run `qiime2 fmt {command-name} --help`


## Examples
Example data can be generated using the `--example-data` flag on each action
described below. This will create a directory structure to match the examples
(you will need to `cd` into the appropriate directory first).

```bash
qiime fmt --example-data fmt-examples/
cd fmt-examples/cc/cc-baseline/
```

