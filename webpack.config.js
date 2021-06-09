const path = require('path');

module.exports = {
  entry: './src/tacoxdna.ts',
  module: {
    rules: [
      {
        test: /\.tsx?$/,
        use: 'ts-loader',
        exclude: /node_modules/,
      },
    ],
  },
  resolve: {
    extensions: ['.tsx', '.ts', '.js'],
  },
  externals: {
    three: {
        commonjs: 'three',
        commonjs2: 'three',
        amd: 'three',
        root: 'THREE',
    },
},
  mode: 'production',
  output: {
    path: path.resolve(__dirname, 'dist'),
    filename: 'tacoxdna.js',
    library: {
        name: 'tacoxdna',
        type: 'umd',
    },
    globalObject: 'this',
  },
};