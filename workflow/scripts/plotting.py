import itertools
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import classification_report, confusion_matrix

plt.ioff()


class ResultsPlotter:
    def __init__(self, y_true, y_pred, target_names, working_dir="."):
        """
        Model characterization functions post-training and testing.

        Args:
        y_true (nparray): 1D npy array containing int values for class
        y_pred (nparray): 1D npy array containing int values for predicted class
        target_names (list[str]): List of string names corresponding to categorical labels 0-n
        working_dir (str/pathlike): Directory to save plots to. Defaults to cwd.
        """
        self.y_true = y_true
        self.y_pred = y_pred
        self.target_names = target_names
        self.working_dir = working_dir

    def plot_confusion_matrix(
        self,
        title="Confusion matrix",
        cmap=None,
        normalize=False,
    ):
        """
        Generates confusion matrix, plots and saves to disk.
        Modified from source: http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html

        Args:
            title (str, optional): Title for plot. Defaults to "Confusion matrix".
            cmap ([type], optional): Colormap. Defaults to None.
            normalize (bool, optional): [description]. Defaults to False.
        """
        cm = self._get_confusion_matrix()

        accuracy = np.trace(cm) / float(np.sum(cm))
        misclass = 1 - accuracy

        if cmap is None:
            cmap = plt.get_cmap("Blues")

        plt.figure(figsize=(8, 8))
        plt.imshow(cm, interpolation="nearest", cmap=cmap)
        plt.title(title)
        plt.colorbar()

        if self.target_names is not None:
            tick_marks = np.arange(len(self.target_names))
            plt.xticks(tick_marks, self.target_names, rotation=45)
            plt.yticks(tick_marks, self.target_names)

        if normalize:
            cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]

        thresh = cm.max() / 1.5 if normalize else cm.max() / 2
        for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
            if normalize:
                plt.text(
                    j,
                    i,
                    "{:0.2f}".format(cm[i, j]),
                    horizontalalignment="center",
                    color="white" if cm[i, j] > thresh else "black",
                )
            else:
                plt.text(
                    j,
                    i,
                    "{:,}".format(cm[i, j]),
                    horizontalalignment="center",
                    color="white" if cm[i, j] > thresh else "black",
                )

        plt.ylabel("True label")
        plt.xlabel(
            "Predicted label\naccuracy={:0.4f}; misclass={:0.4f}".format(
                accuracy, misclass
            )
        )

        plt.savefig(
            os.path.join(
                self.working_dir,
                title.replace(" ", "") + "_discriminator_conf_matrix.png",
            )
        )

    def _get_confusion_matrix(self):
        """Prints confusion matrix, returns for further plotting"""
        conf_mat = confusion_matrix(self.y_true, self.y_pred)

        print("Confusion Matrix")
        print(conf_mat)

        return conf_mat

    def print_classification_report(self):
        """Prints classification report

        Args:
            y_true (nparray): 1D npy array containing int values for class
            y_pred (nparray): 1D npy array containing int values for predicted class
            train_gen (Keras Generator): Training generator used for model training, used for labels
        """
        print("Classification Report")
        print(classification_report(self.y_true, self.y_pred))
