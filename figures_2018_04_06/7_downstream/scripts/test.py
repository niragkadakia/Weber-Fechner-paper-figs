from sklearn.neural_network import MLPClassifier
import sklearn
import inspect
import matplotlib.pyplot as plt
X = [[0., 0., 4, 2, 0, 0, 0, 0], [1., 1.], [3, 5], [0.5, 0.9]]
y = [0, 1, 3, 2]

# wnat: input = 50, hidden = 2000, output = odor #s

clf = MLPClassifier(solver='lbfgs', alpha=1e-5,
                     hidden_layer_sizes=(250), random_state=1)

clf.fit(X, y)                         
print clf.coefs_[1].shape


plt.plot(clf.coefs_[0].T)
plt.show()
plt.pcolormesh(clf.coefs_[0])
plt.colorbar()
plt.show()
plt.pcolormesh(clf.coefs_[1])
plt.colorbar()
plt.show()

vals = clf.predict([[2., 2.], [-1., -2.]])
print vals