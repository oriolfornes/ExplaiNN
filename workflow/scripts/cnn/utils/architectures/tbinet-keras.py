import numpy as np
import h5py
import scipy.io
from sklearn import metrics
import pandas as pd
import os
os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32,gpuarray.preallocate=0.3"
import theano
print(theano.config.device)
from keras.layers import Embedding
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dense, Dropout, Activation, Flatten, Layer, merge, Input, Concatenate, Reshape, concatenate,Lambda,multiply,Permute,Reshape,RepeatVector
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.pooling import GlobalMaxPooling1D
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import optimizers
from keras import backend as K
from keras import regularizers

sequence_input = Input(shape=(35,4))

# Convolutional Layer
output = Conv1D(320,kernel_size=26,padding="valid",activation="relu")(sequence_input)
print(output.shape)
# (?, 975, 320)
output = MaxPooling1D(pool_size=13, strides=13)(output)
output = Dropout(0.2)(output)
exit(0)

#Attention Layer
attention = Dense(1)(output)
attention = Permute((2, 1))(attention)
attention = Activation('softmax')(attention)
attention = Permute((2, 1))(attention)
attention = Lambda(lambda x: K.mean(x, axis=2), name='attention',output_shape=(75,))(attention)
attention = RepeatVector(320)(attention)
attention = Permute((2,1))(attention)
output = multiply([output, attention])

#BiLSTM Layer
output = Bidirectional(LSTM(320,return_sequences=True))(output)
output = Dropout(0.5)(output)

flat_output = Flatten()(output)

#FC Layer
FC_output = Dense(695)(flat_output)
FC_output = Activation('relu')(FC_output)

#Output Layer
output = Dense(690)(FC_output)
output = Activation('sigmoid')(output)

model = Model(inputs=sequence_input, outputs=output)

print('compiling model')
model.compile(loss='binary_crossentropy', optimizer='adam')

print('model summary')
model.summary()
exit(0)
checkpointer = ModelCheckpoint(filepath="./model/tbinet.{epoch:02d}-{val_loss:.2f}.hdf5", verbose=1, save_best_only=False)
earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)

model.fit(X_train, y_train, batch_size=100, epochs=60, shuffle=True, verbose=1, validation_data=(np.transpose(validmat['validxdata'],axes=(0,2,1)),validmat['validdata'][:,125:815]), callbacks=[checkpointer,earlystopper])

model.save('./model/tbinet.h5')