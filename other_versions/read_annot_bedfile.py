from hail.table import Table

def locus_interval_expr(contig, start, end, includes_start, includes_end,
                        reference_genome, skip_invalid_intervals):
    if reference_genome:
        if skip_invalid_intervals:
            is_valid_locus_interval = (
                (hl.is_valid_contig(contig, reference_genome) &
                 (hl.is_valid_locus(contig, start, reference_genome) |
                  (~hl.bool(includes_start) & (start == 0))) &
                 (hl.is_valid_locus(contig, end, reference_genome) |
                  (~hl.bool(includes_end) & hl.is_valid_locus(contig, end - 1, reference_genome)))))

            return hl.or_missing(is_valid_locus_interval,
                                 hl.locus_interval(contig, start, end,
                                                   includes_start, includes_end,
                                                   reference_genome))
        else:
            return hl.locus_interval(contig, start, end, includes_start,
                                     includes_end, reference_genome)
    else:
        return hl.interval(hl.struct(contig=contig, position=start),
                           hl.struct(contig=contig, position=end),
                           includes_start,
                           includes_end)

def import_bed(path, reference_genome='default', skip_invalid_intervals=True) -> Table:
    t = import_table(path, no_header=True, delimiter=r"\s+", impute=False,
                     skip_blank_lines=True, types={'f0': tstr, 'f1': tint32,
                                                   'f2': tint32, 'f3': tstr,
                                                   'f4': tstr},
                     comment=["""^browser.*""", """^track.*""",
                              r"""^\w+=("[\w\d ]+"|\d+).*"""])
    
    t = t.annotate(f0 = hl.str(t.f0).replace('chr',''))
    
    if t.row.dtype == tstruct(f0=tstr, f1=tint32, f2=tint32):
        t = t.select(interval=locus_interval_expr(t['f0'],
                                                  t['f1'] + 1,
                                                  t['f2'],
                                                  True,
                                                  True,
                                                  reference_genome,
                                                  skip_invalid_intervals))

    elif len(t.row) >= 4 and tstruct(**dict([(n, typ) for n, typ in t.row.dtype._field_types.items()][:4])) == tstruct(f0=tstr, f1=tint32, f2=tint32, f3=tstr):
        t = t.select(interval=locus_interval_expr(t['f0'],
                                                  t['f1'] + 1,
                                                  t['f2'],
                                                  True,
                                                  True,
                                                  reference_genome,
                                                  skip_invalid_intervals),
                     target=t['f3'])

    else:
        raise FatalError("too few fields for BED file: expected 3 or more, but found {}".format(len(t.row)))

    if skip_invalid_intervals and reference_genome:
        t = t.filter(hl.is_defined(t.interval))

    return t.key_by('interval')
