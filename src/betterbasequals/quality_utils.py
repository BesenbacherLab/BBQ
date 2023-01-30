from math import sqrt


def get_error_probs(a, b, c):
    """calculate error probabilities for the three different
    mutation types at a type of site (A or C)

    Args:
        a (Float): Rate of incompatibilities of first type fx. n(A->A/C, incomp)/n(A->A/?)
        b (Float): Rate of incompatibilities of second type fx. n(A->A/G, incomp)/n(A->A/?)
        c (Float): Rate of incompatibilities of third type fx. n(A->A/T, incomp)/n(A->A/?)
    """
    # No error prob:
    abc_sum = a + b + c
    d = 1 - 2*abc_sum
    p_AA = 0.5 + sqrt(d)/2
    p_AC = a/(2*p_AA)
    p_AG = b/(2*p_AA)
    p_AT = c/(2*p_AA)
    print(d)
    print(p_AA)
    return p_AC, p_AG, p_AT
sum([x for x in get_error_probs(1/100,1/1000,1/1000)])


def get_error_probs(a, b, c):
    """calculate error probabilities for the three different
    mutation types at a type of site (A or C)

    Args:
        a (Float): Rate of incompatibilities of first type fx. n(A->C, incomp)/n(A->A match)
        b (Float): Rate of incompatibilities of second type fx. n(A->G, incomp)/n(A->A match)
        c (Float): Rate of incompatibilities of third type fx. n(A->T, incomp)/n(A->A match)
    """

    # No error prob:
    abc_sum = a + b + c
    d = 1 - 2*abc_sum
    #print(f'abc_sum={abc_sum} d={d}')
    p_no_error_1 = 1- ((-1 + sqrt(d)) / -2)
    p_no_error_2 = 1- ((-1 - sqrt(d)) / -2)
    #print(p_no_error_1)
    #print(p_no_error_2)
    p_no_error= max(p_no_error_1, p_no_error_2)
    pa = a/(2*p_no_error)
    pb = b/(2*p_no_error)
    pc = c/(2*p_no_error)
    print(pa,pb,pc)
    #print(b/p_no_error)    
    #print(c/p_no_error)
    
