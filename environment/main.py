from code.inference import inference
from code.simulation import simulation
from code.training import training

import click


@click.group()
def cli():
    """A command-line tool with multiple commands."""
    pass


# Add commands to the CLI group
cli.add_command(simulation)
cli.add_command(inference)
cli.add_command(training)

if __name__ == "__main__":
    cli()
