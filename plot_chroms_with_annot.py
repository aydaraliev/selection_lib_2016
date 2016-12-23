for i in range(1,23):
    print("chr %s"%i)
    ijma = ki["%s"%i]
    obiach = ko["%s"%i]
    ijma_peaks = detect_peaks(ijma)
    obiach_peaks = detect_peaks(obiach)
    ijma_leg = plot_stats(ijma, "black", "ijma", peaks_pos = ijma_peaks, genes_to_annot = ki_genes)
    obiach_leg = plot_stats(obiach, "red", "obiach", peaks_pos = obiach_peaks, genes_to_annot = ko_genes)
    annotate_axes("Ijma", "Obiachevo", i, ijma_leg, obiach_leg, "iHS")
    plt.savefig("%s.png"%i)
    plt.clf()
