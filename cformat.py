#!/usr/bin/python
#coding=utf-8
# Created on 
# Author: Zhihua Pei
# Organization: GenePlus
# website: www.zilhua.com
# github: github.com/zilhua
import os
import re
import warnings
sys.dont_write_bytecode = True
sys.path.append(os.path.split(os.path.realpath(__file__))[0])

class VcfInfo(object):
    '''vcf infos: chr pos start end ref alt others'''
    chr, start, end, ref, alt, id, qual, info, strformat, samples, sampleinfos = \
        ['' for i in xrange(11)]
    pass


def vcf2bed(line='', file='', format='4.1', outformat='iterator'):
    _checkvcf2bed(line, file)
    if line:
        return _splitline(line, format)
    else:
        if not os.path.isfile(file):
            raise IOError("input file:{0} not exists!".format(file))
        return _splitfile(file, format, outformat)


def _splitfile(file, format, outformat):
    if outformat != 'iterator':
        out_f = open(outformat, 'w')
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if not line:
                continue
            fs = re.split('\t', line)
            if re.search(r'^##', line):
                continue
            if re.search(r'^#', line):
                samples = '\t'.join(fs[9:])
                continue
            linevcfinfos = _splitline(line, format)
            for l in linevcfinfos:
                l.samples = samples
                if outformat == 'iterator':
                    yield l
                else:
                    out_f.write('\t'.join([l.chr, l.start, l.end, l.ref,
                        l.alt, l.id, l.filter, l.info, l.strformat,
                        l.sampleinfos]) + '\n')
    if outformat != 'iterator':
        out_f.close()


def _splitline(line, format):
    format = _checkformat(line, format)
    if str(format) in ['4.1', '4.0']:
        lines = re.split(os.linesep, line)
        for f in lines:
            fs = re.split('\t', f.strip())
            for alt in re.split(',', fs[4]):
                nstart, nend, nref, nalt = svcf2bed(fs[1], fs[3], alt)
                for i,j,k,m in zip(nstart, nend, nref, nalt):
                    vcfinfos = VcfInfo()
                    vcfinfos.chr, vcfinfos.id, vcfinfos.qual, vcfinfos.filter, \
                        vcfinfos.info, vcfinfos.strformat, vcfinfos.sampleinfos = \
                        fs[0], fs[2], fs[5], fs[6], fs[7], fs[8], '\t'.join(fs[9:])
                    vcfinfos.start, vcfinfos.end, vcfinfos.ref, vcfinfos.alt = \
                        i, j, k, m
                    yield vcfinfos


def svcf2bed(pos, ref, alt):
    '''
    :param pos: vcf pos
    :param ref: vcf ref
    :param alt: vcf alt
    :return: start[list], end[list], ref[list], alt[list]
    '''
    start, end = int(pos), int(pos)
    if len(re.split(',', alt)) > 1:
        raise IOError('not support mutil alt field'.format(alt))
    if len(ref) == len(alt):
        #snv
        start = start -1
        return [start], [end], [ref], [alt]
    elif len(ref) != len(alt) and (len(ref) == 1 or len(alt) ==1):
        #indel
        if len(ref) > len(alt):
            #deletion
            diff = ref[len(alt):]
            return [start], [end+len(diff)], [diff], ['.']
        else:
            #insertion
            diff = alt[len(ref)]
            return [start], [end], ['.'], [diff]
    elif len(ref) != len(alt) and (len(ref) != 1 and len(alt) !=1):
        #delin
        warnings.warn("not suport delins format yet!",SyntaxWarning)
        return [start], [end+len(ref[1:])], [ref[1:]], [alt[1:]]
    else:
        return [], [], [], []


def _checkvcf2bed(line, file):
    if not line and not file:
        raise IOError("No input file or stdin\n")
    if line and os.path.exists(file):
        raise IOError("stdin and file:{0} can't coexists\n".format(file))


def _checkformat(line, format):
    if str(format) in ['4.0', '4.1', '4.2']:
        return str(format)
    raise IOError('{0} format not support yet!'.format(format))


if __name__ == "__main__":
    vcfinfos = vcf2bed(file="170006090BD.combined.vcf", format='4.1', outformat='iterator')
    for vcfinfo in vcfinfos:
        print ">>>>>>>>>start:" + str(vcfinfo.start) + '\tend:' + str(vcfinfo.end)
        print ">>>>>>>>>ref:" + vcfinfo.ref + "\talt:" + vcfinfo.alt
    '''
    with open("170006090BD.combined.vcf", 'r') as f:
        for line in f.readlines():
            if re.search(r'^#', line):
                continue
            vcfinfos = vcf2bed(line,format='4.1')
            for vcfinfo in vcfinfos:
                print ">>>>>>>>>start:" + str(vcfinfo.start)
                print ">>>>>>>>>end:" + str(vcfinfo.end)
                print ">>>>>>>>>ref:" + str(vcfinfo.ref)
                print ">>>>>>>>>alt:" + str(vcfinfo.alt)
                print ">>>>>>>>>infos:" + str(vcfinfo.sampleinfos)
                print ">>>>>>>>>sample:" + str(vcfinfo.samples)
                print ">>>>>>>>>qual:" + str(vcfinfo.qual)
                print ">>>>>>>>>filter:" + str(vcfinfo.filter)
                print ">>>>>>>>>chr:" + str(vcfinfo.chr)
                print ">>>>>>>>>ID:" + str(vcfinfo.id)
                print ">>>>>>>>>ID:" + str(vcfinfo.info)
                print ">>>>>>>>>ID:" + str(vcfinfo.strformat)
                print '\t'.join([vcfinfo.chr, str(vcfinfo.start), str(vcfinfo.end), vcfinfo.ref, vcfinfo.alt, vcfinfo.qual, vcfinfo.sampleinfos, vcfinfo.samples ])
    '''
