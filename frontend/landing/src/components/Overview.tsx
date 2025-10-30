import { Sparkles, Brain, TrendingUp, Zap } from 'lucide-react';

export function Overview() {
  return (
    <div className="max-w-6xl mx-auto space-y-12 animate-in fade-in duration-700">
      {/* Hero Section */}
      <div className="text-center space-y-6">
        <h1 className="text-6xl">
          <span className="bg-gradient-to-r from-purple-400 via-purple-500 to-purple-600 bg-clip-text text-transparent">
            Step into the
          </span>
          <br />
          <span className="text-white">Future</span>
        </h1>
        <p className="text-xl text-purple-300 max-w-3xl mx-auto">
          AI-Driven Drug Discovery Platform
        </p>
      </div>

      {/* Main Overview Card */}
      <div className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-2xl p-10 backdrop-blur-sm shadow-2xl shadow-purple-500/10">
        <div className="flex items-start gap-4 mb-6">
          <div className="p-3 bg-purple-500/20 rounded-lg">
            <Brain className="w-8 h-8 text-purple-400" />
          </div>
          <div>
            <h2 className="text-3xl text-white mb-4">Overview</h2>
            <div className="w-20 h-1 bg-gradient-to-r from-purple-500 to-purple-700 rounded-full"></div>
          </div>
        </div>
        
        <p className="text-lg text-gray-300 leading-relaxed">
          This project is an AI-driven platform designed to help chemists discover and design new drug molecules more efficiently by leveraging machine learning models that predict critical chemical and biological properties. Users upload molecular data and interact with specialized prediction, screening, and ranking pipelines powered by trained models. These tools enable chemists to generate new molecules, screen large candidate libraries, predict key properties such as activity and toxicity, and prioritize compounds for lab testing—dramatically accelerating the discovery process while reducing experimental costs.
        </p>
      </div>

      {/* Key Features Grid */}
      <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
        <div className="bg-gradient-to-br from-purple-900/20 to-purple-950/10 border border-purple-500/20 rounded-xl p-6 hover:border-purple-400/40 transition-all hover:shadow-lg hover:shadow-purple-500/20">
          <div className="p-3 bg-purple-500/10 rounded-lg w-fit mb-4">
            <Sparkles className="w-6 h-6 text-purple-400" />
          </div>
          <h3 className="text-xl text-white mb-3">Generate New Molecules</h3>
          <p className="text-gray-400">
            Create entirely new molecular structures using deep generative AI models
          </p>
        </div>

        <div className="bg-gradient-to-br from-purple-900/20 to-purple-950/10 border border-purple-500/20 rounded-xl p-6 hover:border-purple-400/40 transition-all hover:shadow-lg hover:shadow-purple-500/20">
          <div className="p-3 bg-purple-500/10 rounded-lg w-fit mb-4">
            <TrendingUp className="w-6 h-6 text-purple-400" />
          </div>
          <h3 className="text-xl text-white mb-3">Predict Properties</h3>
          <p className="text-gray-400">
            Estimate activity, toxicity, solubility and drug-likeness with ML models
          </p>
        </div>

        <div className="bg-gradient-to-br from-purple-900/20 to-purple-950/10 border border-purple-500/20 rounded-xl p-6 hover:border-purple-400/40 transition-all hover:shadow-lg hover:shadow-purple-500/20">
          <div className="p-3 bg-purple-500/10 rounded-lg w-fit mb-4">
            <Zap className="w-6 h-6 text-purple-400" />
          </div>
          <h3 className="text-xl text-white mb-3">Accelerate Discovery</h3>
          <p className="text-gray-400">
            Reduce experimental costs while dramatically speeding up the discovery process
          </p>
        </div>
      </div>

      {/* Stats Section */}
      <div className="grid grid-cols-3 gap-8 bg-gradient-to-r from-purple-900/20 via-purple-800/20 to-purple-900/20 border border-purple-500/30 rounded-xl p-8">
        <div className="text-center">
          <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">4</div>
          <div className="text-gray-400">AI Pipelines</div>
        </div>
        <div className="text-center border-l border-r border-purple-500/20">
          <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">∞</div>
          <div className="text-gray-400">Molecules Generated</div>
        </div>
        <div className="text-center">
          <div className="text-4xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">10x</div>
          <div className="text-gray-400">Faster Discovery</div>
        </div>
      </div>
    </div>
  );
}
