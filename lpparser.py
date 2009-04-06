import sys
import tokenizer

class lpparser(object):
    def __init__(self,filename):
        self.tokens = [
            ("variable", "[A-Za-z]+[0-9]+"),
            ("number", "[0-9]+"),
            ("sign", "[\+\-]"),
            ("relational", "[<>]{0,1}="),
            ("eol", "\n"), 
            ]



        self.lex = tokenizer.tokenizer(self.tokens,filename)
        # Prime the machine by reading first token
        self.current = self.lex.next()

    def accept(self,tokenType):
        return self.current[0] == tokenType

    def pgm(self):
        done = False
        while not done:
            self.eol()

    def eol(self):
        self.accept('eol')

    
        
if __name__ == "__main__":
    parser = lpparser("testcube.fm")
    for token in parser.lex.next():
        print token[0],token[1]
    
    parser = lpparser("testcube.fm")
    parser.pgm()
