import matplotlib.pyplot as plt
import numpy as np

def metrics_plot(acc_list, prec_list, recall_list, f1_list, n_folds):
    """
    Define labels and dimensions for subplot displaying metric values for fold subset predictions

    Args:
        acc_list (list): list of model predicition accuracies for each dataset fold after kmeans
        prec_list (list): list of model predicition precision values for each dataset fold after kmeans
        recall_list (list): list of model predicition recall values for each dataset fold after kmeans
        f1_list (list): list of model predicition f1 scores for each dataset fold after kmeans

    Returns:
        avgs (list): list of floats representing mean value of each metric computed for each fold subset
    """

    fig, ax = plt.subplots(figsize=(10,6))

    #define number and types of metrics to plot
    all_metrics = np.array([acc_list, prec_list, recall_list, f1_list])
    num_metrics = len(all_metrics)

    bar_width = 0.15
    indices = np.arange(num_metrics)

    color = plt.cm.viridis(np.linspace(0, 1, all_metrics.shape[1]))

    #define plot dimensions, color, etc.
    for val in range(all_metrics.shape[1]):
        ax.bar(indices + bar_width * val - bar_width * (all_metrics.shape[1] // 2), 
               all_metrics[:, val] * 100, label=f"Fold {val+1}", width = bar_width,
               color = color[val], edgecolor = 'black')

    avgs = np.mean(all_metrics, axis=1)
    
    #define subplot labels and locations
    sub_labels = ['Accuracy', 'Precision', 'Recall', 'f1 Score']
    ax.set_xticks(indices)
    ax.set_xticklabels(sub_labels)

    #set x, y, title
    ax.set_xlabel(f'Metrics across all {n_folds} folds')
    ax.set_ylabel('Performance (%)')
    ax.set_title('Model Performance of Correctly Predicting Cancer Subtype')

    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.show()

    return avgs