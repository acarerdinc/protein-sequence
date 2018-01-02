from keras.models import Model
from keras.layers import Input, LSTM, Dense, Dropout, Masking
from keras.regularizers import l2, l1
import numpy as np
from readdata import read_data
import collections

x_train, y_train, x_text, y_test = read_data(level=0, length_limit=500)
# x_train = x_train * 5
# y_train = y_train * 5

# Vectorize the data.
input_texts = []
target_texts = []
input_characters = set()
targets = set()
for i, (input_text, target) in enumerate(zip(x_train, y_train)):
    input_texts.append(input_text)
    target_texts.append(target)
    for char in input_text:
        if char not in input_characters:
            input_characters.add(char)
    if target not in targets:
        targets.add(target)

input_characters = sorted(list(input_characters))
targets = sorted(list(targets))
num_encoder_tokens = len(input_characters)
nb_targets = len(targets)
max_encoder_seq_length = max([len(txt) for txt in input_texts])
max_decoder_seq_length = max([len(txt) for txt in target_texts])


print('Number of samples:', len(input_texts))
print('Number of unique input tokens:', num_encoder_tokens)
print('Number of unique output classes:', nb_targets)
print('Max sequence length for inputs:', max_encoder_seq_length)
print('Max sequence length for outputs:', max_decoder_seq_length)
print('Class counts:', [(key, value) for
                        key, value in (collections.Counter(y_train)).iteritems()])

input_token_index = dict(
    [(char, i) for i, char in enumerate(input_characters)])
target_token_index = dict(
    [(char, i) for i, char in enumerate(targets)])

# Reverse-lookup token index to decode sequences back to
# something readable.
reverse_input_char_index = dict(
    (i, char) for char, i in input_token_index.items())
reverse_target_index = dict(
    (i, char) for char, i in target_token_index.items())



encoder_input_data = np.zeros(
    (len(input_texts), max_encoder_seq_length, num_encoder_tokens),
    dtype='float32')
y_encoded = np.zeros((len(input_texts), nb_targets), dtype='float32')

for i, (input_text, target_text) in enumerate(zip(input_texts, target_texts)):
    for t, char in enumerate(input_text):
        encoder_input_data[i, -(t+1), input_token_index[char]] = 1.
    y_encoded[i, target_token_index[target_text]] = 1.



batch_size = 32  # Batch size for training.
epochs = 5  # Number of epochs to train for.
latent_dim = 32  # Latent dimensionality of the encoding space.

encoder_inputs = Input(shape=(None, num_encoder_tokens))
mask_1 = Masking(mask_value=0.0)(encoder_inputs)
encoder_1 = LSTM(latent_dim, return_sequences=False,  W_regularizer=l2(0.001))(mask_1)
# encoder_2 = LSTM(latent_dim*2,  W_regularizer=l2(0.001), recurrent_dropout=0.3)(encoder_1)
# encoder_states = [state_h, state_c]
# dense_1 = Dense(64, activation='relu')(encoder_1)
# dropout_1 = Dropout(0.50)(dense_1)
# dense_2 = Dense(128, activation='relu')(dropout_1)
# dropout_2 = Dropout(0.50)(dense_1)
dense_outputs = Dense(nb_targets, activation='softmax')(encoder_1)

model = Model(encoder_inputs, dense_outputs)

# Run training
model.compile(optimizer='rmsprop', loss='categorical_crossentropy', metrics=['accuracy'])
model.summary()
history = model.fit(encoder_input_data, y_encoded,
          batch_size=batch_size,
          epochs=epochs,
          validation_split=0.2,
          verbose=1)

a = 5
print "end"

from keras.utils import CustomObjectScope