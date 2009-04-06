import re

class tokenizer(object):
    def __init__(self, tokens, filename):
        self.tokens = tokens
        self.stream = 0
        parts = []
        for name, part in tokens:
            parts.append("(?P<%s>%s)" % (name, part))
        self.regexpString = "|".join(parts)
        self.regexp = re.compile(self.regexpString, re.MULTILINE)

        try:
            self.stream = open(filename, "r")
        except IOError, e:
            print e

    def cleanup(self):
        self.stream.close()

    def next(self):
        # yield lexemes
        for text in self.stream:
            for match in self.regexp.finditer(text):
                found = False
                for name, rexp in self.tokens:
                    m = match.group(name)
                    if m is not None:
                        yield (name, m)
                        break
        self.stream.close()
        yield ('eof','eof')
