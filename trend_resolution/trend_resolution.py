def init_trend(state):
    if state == 'up':
        # potential trend because utr are sorted by length already 
        pot = -1  # shortening
    elif state == 'down':
        pot = 1 # lengthening
    else:
        raise
    return pot


def define_trend_no_nan(changes):
    """changes should be sorted by utr length, and contain no nan values"""
    trend = 0
    if len(changes) <= 1:
        # there is no point talking about trend if only 1 possibility of 3'UTR length exists
        return trend
    
    transition_count = 0
    state = changes[0]
    trend = init_trend(state)
        
    for i in changes[1:]:
        if i != state:
            transition_count += 1
            state = i
    # print(transition_count)
    # should be exactly 1 transition, if above 1, see CHURC1, LUSC
    if transition_count != 1:
        trend = 0 # reset to undefined trend
    return trend

assert define_trend_no_nan(['up', 'down']) == -1 # shortening
assert define_trend_no_nan(['down', 'up']) == 1  # lengthening
assert define_trend_no_nan(['down']) == 0  # undefined
assert define_trend_no_nan(['up']) == 0  # undefined
assert define_trend_no_nan(['down', 'up', 'down']) == 0  # undefined
assert define_trend_no_nan(['up', 'down', 'up']) == 0  # undefined

def is_nan(val):
    """handles mix of float and string (bad design if you find this function useful)"""
    return True if isinstance(val, float) and np.isnan(val) else False


def define_trend(changes):
    """
    changes should be sorted by utr length, and possibly contain nan values
    
    changes should be partitioned after sorted by utr length, in order to assign a trend,
    otherwise unassigned. e.g.
        down, down, nan, up implies lengthening
        up, down, nan implies lengthening
        
    nan means no significant change, but they have to be there as they're evidence for an 
    alternative 3'UTR with a different length, esp. when only up or down is found for a given stop codon,
    hence for comparison purpose
    """
    if np.unique(changes).shape[0] == 1:
        return 0 # undefined, no trend can be defined
    
    if all([is_nan(_) for _ in changes]):
        # np.unique([np.nan, np.nan]) still has a shape of (2,)
        return 0
    
    if 'up' in changes and 'down' in changes:
        # then ignore nan
        changes = [_ for _ in changes if not is_nan(_)]
    # nan matters when there is only 'down' or 'up', and is 
    # converted to the opposite change, in relatively terms
    elif 'up' in changes:
        changes = [_ if not is_nan(_) else 'down' for _ in changes]
    elif 'down' in changes:
        changes = [_ if not is_nan(_) else 'up' for _ in changes ]
#     print(changes)
    return define_trend_no_nan(changes)

assert define_trend(['up', 'down']) == -1 # shortening
assert define_trend(['down', 'up']) == 1  # lengthening
assert define_trend(['down']) == 0  # undefined
assert define_trend(['down', 'up', 'down']) == 0  # undefined
assert define_trend(['up', 'down', 'up']) == 0  # undefined
# with nan
assert define_trend([np.nan, 'down', 'up']) == 1
assert define_trend(['down', 'down', np.nan]) == 1
assert define_trend(['up', 'down', np.nan]) == -1
assert define_trend(['up', np.nan]) == -1
assert define_trend(['up', 'up', np.nan]) == -1
assert define_trend(['down', np.nan]) == 1
assert define_trend(['down', 'down', np.nan]) == 1
assert define_trend(['down', 'up', np.nan]) == 1
assert define_trend([np.nan, 'up']) == 1
assert define_trend(['up', 'up', np.nan]) == -1
assert define_trend(['down', np.nan]) == 1
assert define_trend(['down', 'down', np.nan]) == 1
assert define_trend([np.nan, np.nan]) == 0

def sc_level_comp(grp_by_sc):
    """stop codon-level comparison"""
    changes = grp_by_sc.sort_values('utr_len').N2T_ratio_change_sig.values.tolist()
#     print(changes)
    return define_trend(changes)

def assign_utr_trend(grp):
    pots = grp.groupby('sc').apply(sc_level_comp)
    return potsdef resolve_trends(trends):
    """resolve possibly multiple trends into one, e.g.
    [shortening, undefined] => shortening
    [lengthening, undefined] => lengthening
    [lengthening, undefined, shortening] => undefined
    [undefined] => undefined
    """
    trends = np.unique(trends)
    if len(trends) == 1:
        return trends[0]
    elif len(trends) == 2:
        if 0 in trends:
            return sum(trends)
        else:
            return 0
    else:
        return 0
        
assert resolve_trends([0]) == 0
assert resolve_trends([1]) == 1
assert resolve_trends([-1]) == -1
assert resolve_trends([-1, 1]) == 0
assert resolve_trends([-1, 0]) == -1
assert resolve_trends([1, 0]) == 1
assert resolve_trends([-1, 0, 1]) == 0


def resolve_trends(trends):
    """resolve possibly multiple trends into one, e.g.
    [shortening, undefined] => shortening
    [lengthening, undefined] => lengthening
    [lengthening, undefined, shortening] => undefined
    [undefined] => undefined
    """
    trends = np.unique(trends)
    if len(trends) == 1:
        return trends[0]
    elif len(trends) == 2:
        if 0 in trends:
            return sum(trends)
        else:
            return 0
    else:
        return 0
        
assert resolve_trends([0]) == 0
assert resolve_trends([1]) == 1
assert resolve_trends([-1]) == -1
assert resolve_trends([-1, 1]) == 0
assert resolve_trends([-1, 0]) == -1
assert resolve_trends([1, 0]) == 1
assert resolve_trends([-1, 0, 1]) == 0
