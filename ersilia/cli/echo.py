try:
    import emoji
except:
    emoji = None
import click


def echo(text, **styles):
    if emoji is not None:
        text = emoji.emojize(text)
    return click.echo(click.style(text, **styles))
