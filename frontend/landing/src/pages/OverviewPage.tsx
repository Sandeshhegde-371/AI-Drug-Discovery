import { Sparkles, Brain, TrendingUp, Zap, Home, ArrowRight, CheckCircle2 } from 'lucide-react';
import { Link } from 'react-router-dom';

export function OverviewPage() {
  return (
    <div className="min-h-screen bg-[#0a0a0f] text-white">
      {/* Gradient Backgrounds */}
      <div className="fixed inset-0 overflow-hidden pointer-events-none">
        <div className="absolute top-0 right-0 w-[600px] h-[600px] bg-purple-600/20 rounded-full blur-[120px]" />
        <div className="absolute bottom-0 left-0 w-[500px] h-[500px] bg-blue-600/15 rounded-full blur-[100px]" />
      </div>

      <div className="relative z-10">
        {/* Header */}
        <header className="px-8 py-6 flex items-center justify-between backdrop-blur-sm bg-white/5 border-b border-white/10 sticky top-0 z-50">
          <Link to="/" className="flex items-center gap-3 hover:opacity-80 transition-opacity group">
            <div className="p-2 bg-purple-500/20 rounded-lg group-hover:bg-purple-500/30 transition-colors">
              <Home className="w-5 h-5 text-purple-400" />
            </div>
            <span className="text-white/80">Back to Home</span>
          </Link>
          
          <div className="flex items-center gap-3">
            <div className="p-2 bg-gradient-to-br from-purple-500 to-purple-700 rounded-lg">
              <Brain className="w-6 h-6 text-white" />
            </div>
            <div className="flex flex-col">
              <span className="text-xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent">
                AI-DRIVEN DDD
              </span>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            <Link to="/uses" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Uses
            </Link>
            <Link to="/sample" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Sample
            </Link>
          </div>
        </header>

        {/* Content */}
        <main className="container mx-auto px-8 py-16">
          <div className="max-w-6xl mx-auto space-y-16">
            {/* Hero Section */}
            <div className="text-center space-y-6">
              <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-purple-500/10 border border-purple-500/20 mb-4">
                <Brain className="w-4 h-4 text-purple-400" />
                <span className="text-sm text-purple-300">Platform Overview</span>
              </div>
              
              <h1 className="text-6xl leading-tight">
                <span className="block text-white mb-2">Step into the</span>
                <span className="bg-gradient-to-r from-purple-400 via-purple-500 to-purple-600 bg-clip-text text-transparent">
                  Future of Discovery
                </span>
              </h1>
              
              <p className="text-xl text-purple-300 max-w-3xl mx-auto">
                AI-Driven Drug Discovery Platform
              </p>
            </div>

            {/* Main Overview Card */}
            <div className="relative bg-gradient-to-br from-purple-950/40 to-purple-900/20 border border-purple-500/30 rounded-3xl p-12 backdrop-blur-sm shadow-2xl overflow-hidden">
              {/* Decorative Elements */}
              <div className="absolute top-0 right-0 w-64 h-64 bg-purple-500/10 rounded-full blur-3xl" />
              <div className="absolute bottom-0 left-0 w-48 h-48 bg-blue-500/10 rounded-full blur-3xl" />
              
              <div className="relative z-10">
                <div className="flex items-start gap-6 mb-8">
                  <div className="p-4 bg-gradient-to-br from-purple-500 to-purple-700 rounded-2xl shadow-lg shadow-purple-500/50">
                    <Brain className="w-10 h-10 text-white" />
                  </div>
                  <div>
                    <h2 className="text-4xl text-white mb-3">Platform Overview</h2>
                    <div className="w-24 h-1.5 bg-gradient-to-r from-purple-500 to-purple-700 rounded-full"></div>
                  </div>
                </div>
                
                <p className="text-lg text-gray-300 leading-relaxed">
                  This project is an AI-driven platform designed to help chemists discover and design new drug molecules more efficiently by leveraging machine learning models that predict critical chemical and biological properties. Users upload molecular data and interact with specialized prediction, screening, and ranking pipelines powered by trained models. These tools enable chemists to generate new molecules, screen large candidate libraries, predict key properties such as activity and toxicity, and prioritize compounds for lab testing—dramatically accelerating the discovery process while reducing experimental costs.
                </p>
              </div>
            </div>

            {/* Key Capabilities */}
            <div>
              <div className="text-center mb-10">
                <h2 className="text-3xl text-white mb-3">Key Capabilities</h2>
                <p className="text-gray-400">Powerful tools to transform your research</p>
              </div>
              
              <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
                <div className="group bg-gradient-to-br from-purple-950/30 to-purple-900/20 border border-purple-500/20 rounded-2xl p-8 hover:border-purple-400/40 transition-all hover:shadow-lg hover:shadow-purple-500/20 hover:-translate-y-1">
                  <div className="p-4 bg-purple-500/10 rounded-xl w-fit mb-6 group-hover:bg-purple-500/20 transition-colors">
                    <Sparkles className="w-8 h-8 text-purple-400" />
                  </div>
                  <h3 className="text-2xl text-white mb-4">Generate New Molecules</h3>
                  <p className="text-gray-400 leading-relaxed">
                    Create entirely new molecular structures using deep generative AI models that explore chemical space beyond human intuition
                  </p>
                </div>

                <div className="group bg-gradient-to-br from-blue-950/30 to-blue-900/20 border border-blue-500/20 rounded-2xl p-8 hover:border-blue-400/40 transition-all hover:shadow-lg hover:shadow-blue-500/20 hover:-translate-y-1">
                  <div className="p-4 bg-blue-500/10 rounded-xl w-fit mb-6 group-hover:bg-blue-500/20 transition-colors">
                    <TrendingUp className="w-8 h-8 text-blue-400" />
                  </div>
                  <h3 className="text-2xl text-white mb-4">Predict Properties</h3>
                  <p className="text-gray-400 leading-relaxed">
                    Estimate activity, toxicity, solubility and drug-likeness with ML models for rapid hypothesis testing
                  </p>
                </div>

                <div className="group bg-gradient-to-br from-indigo-950/30 to-indigo-900/20 border border-indigo-500/20 rounded-2xl p-8 hover:border-indigo-400/40 transition-all hover:shadow-lg hover:shadow-indigo-500/20 hover:-translate-y-1">
                  <div className="p-4 bg-indigo-500/10 rounded-xl w-fit mb-6 group-hover:bg-indigo-500/20 transition-colors">
                    <Zap className="w-8 h-8 text-indigo-400" />
                  </div>
                  <h3 className="text-2xl text-white mb-4">Accelerate Discovery</h3>
                  <p className="text-gray-400 leading-relaxed">
                    Reduce experimental costs while dramatically speeding up the discovery process from months to days
                  </p>
                </div>
              </div>
            </div>

            {/* Benefits */}
            <div className="bg-gradient-to-br from-purple-950/30 to-blue-950/20 border border-purple-500/20 rounded-3xl p-10">
              <h2 className="text-3xl text-white mb-8 text-center">Why Choose Our Platform</h2>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                {[
                  "Access 4 specialized AI pipelines",
                  "Generate unlimited molecular candidates",
                  "Predict critical properties instantly",
                  "Prioritize compounds intelligently",
                  "Reduce lab testing costs by 70%",
                  "Accelerate discovery by 10x"
                ].map((benefit, index) => (
                  <div key={index} className="flex items-center gap-3">
                    <CheckCircle2 className="w-5 h-5 text-purple-400 flex-shrink-0" />
                    <span className="text-gray-300">{benefit}</span>
                  </div>
                ))}
              </div>
            </div>

            {/* Stats Section */}
            <div className="grid grid-cols-3 gap-8 bg-gradient-to-r from-purple-900/30 via-purple-800/30 to-purple-900/30 border border-purple-500/30 rounded-2xl p-10">
              <div className="text-center">
                <div className="text-5xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">4</div>
                <div className="text-gray-400">AI Pipelines</div>
              </div>
              <div className="text-center border-l border-r border-purple-500/20">
                <div className="text-5xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">∞</div>
                <div className="text-gray-400">Molecules Generated</div>
              </div>
              <div className="text-center">
                <div className="text-5xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent mb-2">10x</div>
                <div className="text-gray-400">Faster Discovery</div>
              </div>
            </div>

            {/* CTA */}
            <div className="text-center">
              <Link 
                to="/uses"
                className="inline-flex items-center gap-3 px-8 py-4 rounded-full bg-gradient-to-r from-purple-500 to-purple-700 text-white hover:shadow-lg hover:shadow-purple-500/50 transition-all hover:scale-105"
              >
                <span>Explore Our Pipelines</span>
                <ArrowRight className="w-5 h-5" />
              </Link>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
}
