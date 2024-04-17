using FHist


"""
    chi2, ndf = chi2ndf_ww(h1::Hist1D, h2::Hist1D)

Calculate the chi2 and ndf between two weighted histograms. 
This is useful for comparing two histograms to see if they are statistically compatible. 
    
The chi2 is calculated as ROOT does. For more details see the [ROOT implementation](https://root.cern.ch/doc/master/classTH1.html#ab7d63c7c177ccbf879b5dc31f2311b27).
"""
function chi2ndf_ww(h1::Hist1D, h2::Hist1D)
    if nbins(h1) != nbins(h2)
        error("Histograms have different number of bins")
    end

    chi2 = 0.0
    ndf = nbins(h1) - 1

    bins1 = bincounts(h1)
    bins2 = bincounts(h2)

    errs1 = sumw2(h1)
    errs2 = sumw2(h2)

    sum1 = sum(bins1)
    sum2 = sum(bins2)

    for (cnt1, cnt2, e1sq, e2sq) in zip(bins1, bins2, errs1, errs2)
        if cnt1 * cnt1 == 0.0 && cnt2 * cnt2 == 0.0
            ndf -= 1
            continue
        end

        if e1sq == 0.0 && e2sq == 0.0
            error("Both errors are zero")
        end

        sigma = sum1 * sum1 * e2sq + sum2 * sum2 * e1sq
        delta = sum2 * cnt1 - sum1 * cnt2
        chi2 += delta * delta / sigma
    end

    return chi2, ndf
end

function chi2test(h1::Hist1D{T}, h2::Hist1D{T}, options::Symbol...) where {T<:Real}


    nbin1 = nbins(h1)
    nbin2 = nbins(h2)

    if nbin1 != nbin2
        error("Histograms have different number of bins")
    end

    chi2 = 0.0
    ndf = nbin1 - 1

    scaled_histogram = :NORM in options
    comparison_uu = :UU in options
    comparison_ww = :WW in options
    comparison_uw = :UW in options

    if scaled_histogram && !comparison_uu
        error("NORM option should be used together with UU option")
    end

    sum1 = 0
    sum2 = 0

    _sumw1 = 0
    _sumw2 = 0

    if comparison_uu && scaled_histogram
        cnt1s = bincounts(h1)
        cnt2s = bincounts(h2)

        _sumw1 = sumw2(h1)
        _sumw2 = sumw2(h2)

        for i in 1:nbin1
            if _sumw1[i] > 0
                cnt1s[i] = floor(cnt1s[i] * cnt1s[i] / _sumw1[i] + 0.5)
            else
                cnt1s[i] = 0.0
            end

            if _sumw2[i] > 0
                cnt2s[i] = floor(cnt2s[i] * cnt2s[i] / _sumw2[i] + 0.5)
            else
                cnt2s[i] = 0.0
            end
        end

        sum1 = sum(cnt1s)
        sum2 = sum(cnt2s)

        if sumw1 <= 0 || sumw2 <= 0
            error("Cannot use option NORM when one histogram has all zero errors")
        end
    else
        sum1 = sum(bincounts(h1))
        sum2 = sum(bincounts(h2))

        if comparison_ww
            _sumw1 = sumw2(h1)
        end

        if comparison_ww || comparison_uw
            _sumw2 = sumw2(h2)
        end
    end

    if comparison_uu
        error("Not implemented yet")
    elseif comparison_uw
        error("Not implemented yet")
    elseif comparison_ww

        bins1 = bincounts(h1)
        bins2 = bincounts(h2)

        errs1 = sumw2(h1)
        errs2 = sumw2(h2)

        for (cnt1, cnt2, e1sq, e2sq) in zip(bins1, bins2, errs1, errs2)
            if cnt1 * cnt1 == 0.0 && cnt2 * cnt2 == 0.0
                ndf -= 1
                continue
            end

            if e1sq == 0.0 && e2sq == 0.0
                error("Both errors are zero")
            end

            sigma = sum1 * sum1 * e2sq + sum2 * sum2 * e1sq
            delta = sum2 * cnt1 - sum1 * cnt2
            chi2 += delta * delta / sigma
        end
    end

    return chi2, ndf
end