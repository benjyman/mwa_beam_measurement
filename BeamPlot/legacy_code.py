     
def readht(path, line=-1, nchars=0, offset=0, chopnl=True, retnull=''):
    """Reads from line (head>=0, tail<0) nchars with offset, all \n cut if chopnl,
    errors or insufficient chars/lines return '' if retnull else error raised"""
    request, find = ('tail',line) if line<0 else ('head',-1*line+1)
    p = Popen([request, str(line), path],stdout=PIPE,stderr=PIPE).communicate()
    out = p[0].strip().split('\n') if chopnl else p[0]
    if p[1] or len(out)<-1*find or len(out[line])<nchars+offset:
        if retnull is False:
           err('Error: '+('%s'%p[1] if p[1] else 'fewer lines/chars than req.'))
        return retnull
    return out[line][offset:offset+nchars] if nchars else out[line]


def csvparse(path, form=None, req=None, head=True, delim=','):
    """Generator of header/line dicts (list header if absent.) output as such or
    as defined by form function if conditions in req function are satisfied"""
    csv    = (line.strip().split(delim) for line in open(path))
    head   = csv.next() if head is True else head
    data   = [d for d in (dict(zip(head,line)) for line in csv) if req is None or req(d)]
    return data if form is None else (form(d) for d in data)



def err(description):
    """Raises attribute error if characteristics are invalid"""
    print os.getcwd
    raise AttributeError(description)



#######################################################
"""
days = [('Tue','24'),('Wed','25'),('Thu','26'),('Fri','27'),('Sat','28')]

def clean(filepath,forcedate=0):
    newfile=ren(filepath,forcedate)
    if 'log' in filepath:
        os.system('cp %s ../Obs/Check/%s'%(filepath,newfile))
        return
    i,errors=0,0
    for line in open(f):
        i+=1
        if not line.strip():
            continue
        row=line.strip()
        s=row.split('$Sp')
        try:
            t, d = float(s[0]), [ord(s[1][j]) for j in range(len(s[1]))]
            if len(row) in [127,128] and len(s[1])==112:
                open('../Obs/Check/%s'%newfile,'a').write('%.2f'%t+s[1]+'\n')
                continue
        except (ValueError,TypeError,IndexError):
            pass
        print '%d:%s'%(i,line.strip())
        errors+=1
    print i,errors
    return

def ren(name,forcedate):
        '''TTTRRXX-20151104-136800138200-NOV15.txt'''
        
        name=name.split('/')[-1]
        string=''
        if 'log' in name:
            string+='log-'
        elif not 'RFEdata' in name:
            return
        if 'rec' in name:
            string+='0'+name.split('-')[-1][0:2]
            string+='YY'+name.split('rec')[1][0:2]+'-'
        elif 'ref' in name:
            n=name.rstrip('.txt')[-1]
            string+='rf'+n+'XX00-'
        if forcedate:
            string+=forcedate+'-'
        else:
            for u,v in days:
                if u in name:
                    string+='201511'+v+'-'
        string+='01368000138200-NOV15.txt'
        return string

initialpath='/media/monolith/SSD240/CERBERUS/'
for f in glob(initialpath + '*'):
    clean(f,forcedate='')

def chform():
    for f in glob('./TimeRows/*'):
	filename=f.rstrip('trows.txt')+'frows.txt'
	filename=filename.lstrip('./TimeRows/')
	times=[]
	g=open(f)
	for l in g:
	    times.append(l[:13])
	g.seek(0)
	h=open(filename,'a')
	for t in range(len(times)):
	    h.write(times[t])
	    if t!=len(times)-1:
		h.write(',')
	    else:
		h.write('\n')
	for i in range(112):
	    lin,chan=len(times),''
	    g.seek(0)
	    for l in g:
		chan+=l[13+i]
	    h.write(chan+'\n')
	g.close()
	h.close()
"""
