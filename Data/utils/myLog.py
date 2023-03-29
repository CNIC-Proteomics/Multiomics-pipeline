from datetime import datetime

def myLog(path):
    
    def logw(msg):
        with open(path, 'a') as f:
            print(msg)
            f.write(f'{datetime.now()} - {msg}\n')
    
    return logw