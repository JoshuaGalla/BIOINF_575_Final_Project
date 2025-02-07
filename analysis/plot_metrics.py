import matplotlib.pyplot as plt
import numpy as np

def metrics_plot(acc_list, prec_list, recall_list, f1_list, auc_list, n_folds):
    """
    """

    fig, ax = plt.subplots(figsize=(10,6))

    all_metrics = np.array([acc_list, prec_list, recall_list, f1_list, auc_list])
    num_metrics = len(all_metrics)

    bar_width = 0.15
    indices = np.arange(num_metrics)

    color = plt.cm.viridis(np.linspace(0, 1, all_metrics.shape[1]))

    for val in range(num_metrics):
        ax.bar(indices + bar_width * val - bar_width * (all_metrics.shape[1] // 2), 
               all_metrics[:, val] * 100, label=f"Fold {val+1}", width = bar_width,
               color = color[val], edgecolor = 'black')

    avgs = np.mean(all_metrics, axis=1)
    
    sub_labels = ['Accuracy', 'Precision', 'Recall', 'f1 Score', 'ROC_AUC Score']
    ax.set_xticks(indices)
    ax.set_xticklabels(sub_labels)

    ax.set_xlabel(f'Metrics across all {n_folds} folds')
    ax.set_ylabel('Performance (%)')
    ax.set_title('Model Performance of Correctly Predicting Cancer Subtype')

    ax.legend()

    plt.tight_layout()
    plt.show()

    return avgs